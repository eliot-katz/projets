#include <rw/pricing/mc_pricer.hpp>

#include <chrono>    // steady_clock, duration_cast
#include <cmath>     // std::exp, std::sqrt
#include <algorithm> // std::min
#include <rw/pricing/analytic_bs.hpp> // prix BS pour E[Y]

namespace rw {
namespace pricing {

namespace {
inline double half_width_from_se(double se) {
  // Z_95 bilatéral
  static constexpr double Z95 = 1.959963984540054;
  return Z95 * se;
}
} // anonymous

MonteCarloResult McPricer::price_european(const rw::market::MarketData& mkt,
                                          const rw::market::Instrument& inst) const {
  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();

  const double T = inst.T;
  const double K = inst.K;

  // Court-circuit T == 0 : payoff déterministe, pas d'actualisation.
  if (T == 0.0) {
    const double ST = mkt.S0;
    const bool is_call = (inst.type == rw::market::OptionType::Call);
    const double payoff = is_call
        ? rw::payoffs::payoff_call(ST, K)
        : rw::payoffs::payoff_put (ST, K);

    MonteCarloResult res{};
    res.price        = payoff;
    res.std_error    = 0.0;
    res.ci_low       = payoff;
    res.ci_high      = payoff;
    res.n_effective  = 0;
    res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - t0).count();
    return res;
  }

  // Préparation RNG, constantes de modèle
  rw::core::NormalRng rng(cfg_.seed);
  const auto& p = model_.params(); // r, q, sigma
  const double r = p.r, q = p.q, sigma = p.sigma;

  const double df = std::exp(-r * T);

  const std::size_t nSteps = (cfg_.n_steps > 0 ? cfg_.n_steps : 1);
  const double dt = (nSteps == 1) ? T : (T / static_cast<double>(nSteps));

  // Helpers payoff
  auto payoff_from_ST = [&](double ST, const rw::market::Instrument& ins) -> double {
    const bool is_call = (ins.type == rw::market::OptionType::Call);
    return is_call ? rw::payoffs::payoff_call(ST, ins.K)
                   : rw::payoffs::payoff_put (ST, ins.K);
  };

  // --------------------- E[Y] analytique (prix BS du CV) ---------------------
  rw::market::Instrument cv_ins = (cfg_.cv_mode == rw::config::CvMode::SameInstrument)
                                    ? inst
                                    : cfg_.cv_instrument;

  // Par prudence, si CustomInstrument avec T différent, on utilise le T du CV tel quel.
  const double EY = (cv_ins.type == rw::market::OptionType::Call)
          ? rw::pricing::price_call_bs(mkt.S0, cv_ins.K, r, q, sigma, cv_ins.T)
          : rw::pricing::price_put_bs (mkt.S0, cv_ins.K, r, q, sigma, cv_ins.T);

  // --------------------- Phase pilote pour β (si demandée) -------------------
  double beta = 0.0;
  std::size_t paths_consumed_pilot = 0; // chemins physiques consommés par le pilote
  if (cfg_.use_control_variate && cfg_.cv_beta_mode == rw::config::BetaMode::Pilot) {
  // Budget pilote en chemins physiques
  std::size_t budget_pilot = std::min(cfg_.cv_pilot_paths, cfg_.n_paths_target);
  if (cfg_.use_antithetic && (budget_pilot % 2u == 1u)) {
    --budget_pilot; // pair pour des paires complètes
  }

  // Pré-calc des incréments
  const double mu_T   = (r - q - 0.5 * sigma * sigma) * T;
  const double vol_T  = sigma * std::sqrt(T);
  const double mu_dt  = (r - q - 0.5 * sigma * sigma) * dt;
  const double vol_dt = sigma * std::sqrt(dt);

  std::size_t n = 0; // nombre d'unités statistiques (paires si anti, chemins sinon)

  if (cfg_.use_antithetic) {
    // *** Nouvel estimateur "compatible anti" ***
    // On estime β sur les SOMMES de paires : SX = X+ + X-, SY = Y+ + Y-
    long double sumSX = 0.0L, sumSY = 0.0L, sumSX_SY = 0.0L, sumSY2 = 0.0L;
    std::size_t consumed = 0;

    while (consumed + 2 <= budget_pilot) {
      double STp, STm;

      if (nSteps == 1) {
        const double Z = rng.sample();
        STp = mkt.S0 * std::exp(mu_T + vol_T * Z);
        STm = mkt.S0 * std::exp(mu_T - vol_T * Z);
      } else {
        double Sp = mkt.S0, Sm = mkt.S0;
        for (std::size_t k = 0; k < nSteps; ++k) {
          const double Zk = rng.sample();
          Sp *= std::exp(mu_dt + vol_dt * Zk);
          Sm *= std::exp(mu_dt - vol_dt * Zk);
        }
        STp = Sp; STm = Sm;
      }

      const double Xp = df * payoff_from_ST(STp, inst);
      const double Xm = df * payoff_from_ST(STm, inst);
      const double Yp = df * payoff_from_ST(STp, cv_ins);
      const double Ym = df * payoff_from_ST(STm, cv_ins);

      const long double SX = static_cast<long double>(Xp + Xm);
      const long double SY = static_cast<long double>(Yp + Ym);

      sumSX    += SX;
      sumSY    += SY;
      sumSX_SY += SX * SY;
      sumSY2   += SY * SY;

      consumed += 2;
      ++n; // une PAIRE = 1 unité
    }

    paths_consumed_pilot = consumed;

    if (n >= 2) {
      const long double EX  = sumSX / n;
      const long double EYp = sumSY / n;
      const long double EXY = sumSX_SY / n;
      const long double EY2 = sumSY2 / n;

      const long double covXY = EXY - EX * EYp;
      const long double varY  = EY2 - EYp * EYp;

      beta = (varY > 0.0L) ? static_cast<double>(covXY / varY) : 0.0;
    } else {
      beta = 0.0;
    }

  } else {
    // *** Cas non-antithétique : estimateur classique sur (X, Y) ***
    long double sumX = 0.0L, sumY = 0.0L, sumXY = 0.0L, sumY2 = 0.0L;
    std::size_t consumed = 0;

    while (consumed + 1 <= budget_pilot) {
      double ST;
      if (nSteps == 1) {
        const double Z = rng.sample();
        ST = mkt.S0 * std::exp(mu_T + vol_T * Z);
      } else {
        double S = mkt.S0;
        for (std::size_t k = 0; k < nSteps; ++k) {
          const double Zk = rng.sample();
          S *= std::exp(mu_dt + vol_dt * Zk);
        }
        ST = S;
      }

      const double X = df * payoff_from_ST(ST, inst);
      const double Y = df * payoff_from_ST(ST, cv_ins);

      sumX  += X;
      sumY  += Y;
      sumXY += static_cast<long double>(X) * static_cast<long double>(Y);
      sumY2 += static_cast<long double>(Y) * static_cast<long double>(Y);

      consumed += 1;
      ++n;
    }

    paths_consumed_pilot = consumed;

    if (n >= 2) {
      const long double Xbar = sumX / n;
      const long double Ybar = sumY / n;
      const long double EXY  = sumXY / n;
      const long double EY2  = sumY2 / n;

      const long double covXY = EXY - Xbar * Ybar;
      const long double varY  = EY2 - Ybar * Ybar;

      beta = (varY > 0.0L) ? static_cast<double>(covXY / varY) : 0.0;
    } else {
      beta = 0.0;
    }
    }
  } else if (cfg_.use_control_variate && cfg_.cv_beta_mode == rw::config::BetaMode::Fixed) {
    beta = cfg_.cv_beta_fixed;
  } else {
    beta = 0.0; // pas de CV
  }

  // --------------------- Phase principale (avec/ sans CV, Anti Ok) ----------
  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) {
    log.reserve((cfg_.n_paths_target + (cfg_.batch_size ? cfg_.batch_size : 1) - 1) /
                (cfg_.batch_size ? cfg_.batch_size : 1));
  }

  auto accumulate_log = [&](std::size_t n_cum_paths) {
    if (!log_enabled_) return;
    const double se = acc.std_error();
    const double half = half_width_from_se(se);
    log.push_back(ConvergencePoint{ n_cum_paths, acc.mean(), half });
  };

  // Reste de budget (chemins physiques) après le pilote
  std::size_t Ncum_paths = paths_consumed_pilot;
  std::size_t budget_remaining = (cfg_.n_paths_target > Ncum_paths) ? (cfg_.n_paths_target - Ncum_paths) : 0;

  const double mu_T   = (r - q - 0.5 * sigma * sigma) * T;
  const double vol_T  = sigma * std::sqrt(T);
  const double mu_dt  = (r - q - 0.5 * sigma * sigma) * dt;
  const double vol_dt = sigma * std::sqrt(dt);

  if (cfg_.use_antithetic) {
    while (budget_remaining >= 2) {
      const std::size_t max_pairs_by_budget = budget_remaining / 2;
      const std::size_t max_pairs_by_batch  = (cfg_.batch_size >= 2) ? (cfg_.batch_size / 2) : max_pairs_by_budget;
      const std::size_t pairs = std::min<std::size_t>(max_pairs_by_budget, max_pairs_by_batch);
      if (pairs == 0) break;

      for (std::size_t i = 0; i < pairs; ++i) {
        double STp, STm;

        if (nSteps == 1) {
          const double Z = rng.sample();
          STp = mkt.S0 * std::exp(mu_T + vol_T * Z);
          STm = mkt.S0 * std::exp(mu_T - vol_T * Z);
        } else {
          double Sp = mkt.S0, Sm = mkt.S0;
          for (std::size_t k = 0; k < nSteps; ++k) {
            const double Zk = rng.sample();
            Sp *= std::exp(mu_dt + vol_dt * Zk);
            Sm *= std::exp(mu_dt - vol_dt * Zk);
          }
          STp = Sp; STm = Sm;
        }

        const double Xp = df * payoff_from_ST(STp, inst);
        const double Xm = df * payoff_from_ST(STm, inst);

        double sample_value;
        if (cfg_.use_control_variate) {
          const double Yp = df * payoff_from_ST(STp, cv_ins);
          const double Ym = df * payoff_from_ST(STm, cv_ins);
          const double X_pair = 0.5 * (Xp + Xm);
          const double Y_pair = 0.5 * (Yp + Ym);
          sample_value = X_pair - beta * (Y_pair - EY);
        } else {
          sample_value = 0.5 * (Xp + Xm); // antithétique simple
        }

        acc.add(sample_value);
      }

      Ncum_paths      += 2 * pairs;
      budget_remaining = (cfg_.n_paths_target > Ncum_paths) ? (cfg_.n_paths_target - Ncum_paths) : 0;

      accumulate_log(Ncum_paths);

      const double se_now = acc.std_error();
      if (cfg_.tolerance > 0.0 && se_now < cfg_.tolerance) break;
    }
  } else {
    while (budget_remaining >= 1) {
      std::size_t batchN = cfg_.batch_size ? std::min(cfg_.batch_size, budget_remaining) : budget_remaining;
      if (batchN == 0) break;

      for (std::size_t i = 0; i < batchN; ++i) {
        double ST;

        if (nSteps == 1) {
          const double Z = rng.sample();
          ST = mkt.S0 * std::exp(mu_T + vol_T * Z);
        } else {
          double S = mkt.S0;
          for (std::size_t k = 0; k < nSteps; ++k) {
            const double Zk = rng.sample();
            S *= std::exp(mu_dt + vol_dt * Zk);
          }
          ST = S;
        }

        const double X = df * payoff_from_ST(ST, inst);
        if (cfg_.use_control_variate) {
          const double Y = df * payoff_from_ST(ST, cv_ins);
          const double Xstar = X - beta * (Y - EY);
          acc.add(Xstar);
        } else {
          acc.add(X);
        }
      }

      Ncum_paths      += batchN;
      budget_remaining = (cfg_.n_paths_target > Ncum_paths) ? (cfg_.n_paths_target - Ncum_paths) : 0;

      accumulate_log(Ncum_paths);

      const double se_now = acc.std_error();
      if (cfg_.tolerance > 0.0 && se_now < cfg_.tolerance) break;
    }
  }

  const auto t1 = clock::now();

  // Résultats finaux
  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  const auto ci    = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = acc.count(); // nb d'échantillons (paires si anti)
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

} // namespace pricing
} // namespace rw
