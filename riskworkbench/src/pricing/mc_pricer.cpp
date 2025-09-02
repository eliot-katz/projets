#include <rw/pricing/mc_pricer.hpp>

#include <chrono> // steady_clock, duration_cast
#include <cmath>  // std::exp

namespace rw {
namespace pricing {

MonteCarloResult McPricer::price_european(const rw::market::MarketData& mkt,
                                          const rw::market::Instrument& inst) const {
  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();

  // Raccourci T == 0 : payoff déterministe, pas d'aléa ni d'actualisation.
  const double T = inst.T;
  const double K = inst.K;
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
    // Journal vide (pas de simulation)
    return res;
  }

  // Préparation RNG et constantes
  rw::core::NormalRng rng(cfg_.seed);
  const auto& p = model_.params(); // r, q, sigma (sigma utilisé par le modèle)
  const double df = std::exp(-p.r * T);

  const std::size_t nSteps = (cfg_.n_steps > 0 ? cfg_.n_steps : 1);
  const double dt = (nSteps == 1) ? T : (T / static_cast<double>(nSteps));

  // Accumulateur streaming
  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) {
    log.reserve((cfg_.n_paths_target + (cfg_.batch_size ? cfg_.batch_size : 1) - 1) /
                (cfg_.batch_size ? cfg_.batch_size : 1));
  }

  std::size_t Ncum = 0;
  while (true) {
    // Taille du lot courant
    std::size_t remaining = (cfg_.n_paths_target > Ncum) ? (cfg_.n_paths_target - Ncum) : 0;
    if (remaining == 0) break;
    std::size_t batchN = cfg_.batch_size;
    if (batchN == 0) batchN = remaining;             // robustesse si batch_size==0
    if (batchN > remaining) batchN = remaining;      // borner le dernier lot

    // Simulation du lot
    for (std::size_t i = 0; i < batchN; ++i) {
      double ST;
      if (nSteps == 1) {
        ST = model_.sample_ST(mkt.S0, T, rng);
      } else {
        double S = mkt.S0;
        for (std::size_t k = 0; k < nSteps; ++k) {
          model_.step_inplace(S, dt, rng);
        }
        ST = S;
      }

      const bool is_call = (inst.type == rw::market::OptionType::Call);
      const double payoff = is_call
          ? rw::payoffs::payoff_call(ST, K)
          : rw::payoffs::payoff_put (ST, K);

      const double pv = df * payoff;
      acc.add(pv);
    }

    Ncum += batchN;

    // Journal de convergence (par lot)
    if (log_enabled_) {
      const double se = acc.std_error();
      const auto ci = rw::core::confidence_interval_95(acc.mean(), se, acc.count());
      const double half = 0.5 * (ci.high - ci.low);
      log.push_back(ConvergencePoint{ Ncum, acc.mean(), half });
    }

    // Critère d'arrêt (erreur standard)
    const double se_now = acc.std_error(); // NaN si n<2 => comparaison false
    if (cfg_.tolerance > 0.0 && se_now < cfg_.tolerance) {
      break;
    }
    if (Ncum >= cfg_.n_paths_target) {
      break;
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
  res.n_effective  = acc.count();
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

} // namespace pricing
} // namespace rw
