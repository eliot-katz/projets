#include "McWorker.hpp"

#include <rw/core/stats.hpp>
#include <rw/payoffs/vanilla.hpp>
#include <rw/pricing/greeks.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

#include <QDebug>
#include <QString>
using rw::core::RunningStats;

namespace {
// Helpers locaux (pas dans le header)
inline double discount(double r, double T) noexcept { return std::exp(-r * T); }

inline double payoff_from_ST(double ST, const rw::market::Instrument& inst) noexcept {
  using rw::market::OptionType;
  return (inst.type == OptionType::Call)
           ? rw::payoffs::payoff_call(ST, inst.K)
           : rw::payoffs::payoff_put (ST, inst.K);
}
} // namespace

namespace gui {

McWorker::McWorker(QObject* parent) : QObject(parent) {}

void McWorker::requestStop() { stop_.store(true, std::memory_order_relaxed); }

void McWorker::runPricing(rw::market::MarketData mkt,
                          rw::market::Instrument inst,
                          rw::config::McConfig cfg,
                          double sigma)
{
  stop_.store(false, std::memory_order_relaxed);

  try {
    // Guards
    if (mkt.S0 <= 0.0) { emit failed("Invalid S0"); return; }
    if (inst.K <= 0.0 || inst.T < 0.0) { emit failed("Invalid K/T"); return; }
    if (sigma < 0.0) { emit failed("Invalid sigma"); return; }

    // Modèle
    const double sig = effectiveSigma_(mkt, inst, sigma, /*applyVega=*/false);
    rw::models::Gbm model({mkt.r, mkt.q, sig});

    qDebug() << "[runPricing]"
         << "S0=" << mkt.S0 << "r=" << mkt.r << "q=" << mkt.q
         << "K=" << inst.K << "T=" << inst.T
         << "sigmaGlobal=" << sigma
         << "smileEnabled=" << smileEnabled_ << "vega%=" << vegaStressPct_;

    qDebug() << "[runPricing]" << "sigmaEff=" << sig;

    const double T  = inst.T;
    const double df = discount(mkt.r, T);

    const auto t0 = std::chrono::steady_clock::now();

    RunningStats acc;
    std::size_t Ncum = 0;
    rw::core::NormalRng rng(cfg.seed);

    const std::size_t nSteps = (cfg.n_steps == 0 ? 1 : cfg.n_steps);
    const double dt = (nSteps <= 1 ? T : (T / static_cast<double>(nSteps)));

    // Cas T==0 : déterministe (pas de discount non plus, car T=0)
    if (T == 0.0) {
      const double ST = mkt.S0;
      const double pv = payoff_from_ST(ST, inst);
      acc.add(pv);
      Ncum = 1;
      const double se   = acc.std_error(); // 0
      const double half = 1.959963984540054 * se;
      emit progress(Ncum, acc.mean(), half);
      const auto t1 = std::chrono::steady_clock::now();
      emit finished(acc.mean(), se, acc.mean()-half, acc.mean()+half, acc.count(),
                    std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
      return;
    }

    // Boucle par lots
    while (Ncum < cfg.n_paths_target && !stop_.load(std::memory_order_relaxed)) {
      const std::size_t budget = cfg.n_paths_target - Ncum;
      const std::size_t batchN = std::min(cfg.batch_size, budget);

      for (std::size_t i = 0; i < batchN; ++i) {
        double ST = 0.0;
        if (nSteps <= 1) {
          ST = model.sample_ST(mkt.S0, T, rng);
        } else {
          double S = mkt.S0;
          for (std::size_t k = 0; k < nSteps; ++k) model.step_inplace(S, dt, rng);
          ST = S;
        }
        const double pv = df * payoff_from_ST(ST, inst);
        acc.add(pv);
      }

      Ncum += batchN;

      // Progress live
      const double se   = acc.std_error();
      const double half = 1.959963984540054 * se;
      emit progress(Ncum, acc.mean(), half);

      // Arrêt sur tolérance (SE)
      if (cfg.tolerance > 0.0 && se < cfg.tolerance) break;
    }

    const auto t1 = std::chrono::steady_clock::now();

    if (stop_.load(std::memory_order_relaxed)) {
      emit canceled();
      return;
    }

    const double se   = acc.std_error();
    const double half = 1.959963984540054 * se;
    emit finished(acc.mean(), se, acc.mean() - half, acc.mean() + half, acc.count(),
                  std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
  } catch (const std::exception& e) {
    emit failed(QString::fromUtf8(e.what()));
  } catch (...) {
    emit failed("Unknown error in McWorker::runPricing");
  }
}

void McWorker::runGreeks(rw::market::MarketData mkt,
                         rw::market::Instrument inst,
                         rw::config::McConfig mc,
                         rw::config::GreeksConfig gk,
                         double sigma,
                         QString method)
{
  stop_.store(false, std::memory_order_relaxed);

  try {
    if (mkt.S0 <= 0.0) { emit failed("Invalid S0"); return; }
    if (inst.K <= 0.0 || inst.T < 0.0) { emit failed("Invalid K/T"); return; }
    if (sigma <= 0.0) { emit failed("Invalid sigma for Greeks (<=0)"); return; }

    qDebug() << "[runGreeks]"
         << "S0=" << mkt.S0 << "r=" << mkt.r << "q=" << mkt.q
         << "K=" << inst.K << "T=" << inst.T
         << "sigmaGlobal=" << sigma
         << "smileEnabled=" << smileEnabled_ << "vega%=" << vegaStressPct_;
    const double sig = effectiveSigma_(mkt, inst, sigma, /*applyVega=*/false);
    qDebug() << "[runGreeks]" << "sigmaEff=" << sig;

    // Modèle & estimator
    rw::models::Gbm model({mkt.r, mkt.q, sig});
    // gk vient de l’UI (bump, anti, crn, steps…)
    rw::pricing::GreeksEstimator est(model, gk, mc);

    const auto t0 = std::chrono::steady_clock::now();

    // On calcule toutes les lettres selon la "méthode" voulue
    double Delta = std::numeric_limits<double>::quiet_NaN();
    double Vega  = std::numeric_limits<double>::quiet_NaN();
    double Gamma = std::numeric_limits<double>::quiet_NaN();
    double Rho   = std::numeric_limits<double>::quiet_NaN();
    double Theta = std::numeric_limits<double>::quiet_NaN();

    const QString m = method.trimmed().toLower();

    if (m == "brv") {
      Delta = est.delta_brv(mkt, inst).price;
      Vega  = est.vega_brv (mkt, inst).price;
      Rho   = est.rho_brv  (mkt, inst).price;
      Theta = est.theta_brv(mkt, inst).price;

      // Gamma : on peut l’obtenir via ΔPW central même en mode BRV (plus stable)
      rw::config::GreeksConfig gk_pw = gk; gk_pw.use_pathwise = true;
      rw::pricing::GreeksEstimator est_pw(model, gk_pw, mc);
      Gamma = est_pw.gamma_from_delta_pw(mkt, inst).price;

    } else if (m == "pw") {
      // Delta pathwise + Gamma via ΔPW central
      Delta = est.delta_pathwise(mkt, inst).price;
      Gamma = est.gamma_from_delta_pw(mkt, inst).price;

      // Le reste en BRV
      Vega  = est.vega_brv (mkt, inst).price;
      Rho   = est.rho_brv  (mkt, inst).price;
      Theta = est.theta_brv(mkt, inst).price;

    } else if (m == "lrm") {
      // Vega via LRM (sans bump) + autres en BRV
      Vega  = est.vega_lrm (mkt, inst).price;
      Delta = est.delta_brv(mkt, inst).price;
      Rho   = est.rho_brv  (mkt, inst).price;
      Theta = est.theta_brv(mkt, inst).price;

      // Gamma via ΔPW central (robuste)
      rw::config::GreeksConfig gk_pw = gk; gk_pw.use_pathwise = true;
      rw::pricing::GreeksEstimator est_pw(model, gk_pw, mc);
      Gamma = est_pw.gamma_from_delta_pw(mkt, inst).price;

    } else {
      emit failed("Unknown greeks method (use brv|pw|lrm)");
      return;
    }

    const auto t1 = std::chrono::steady_clock::now();
    const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    if (stop_.load(std::memory_order_relaxed)) {
      emit canceled();
      return;
    }
    emit greeksFinished(Delta, Vega, Gamma, Rho, Theta, ms);

  } catch (const std::exception& e) {
    emit failed(QString::fromUtf8(e.what()));
  } catch (...) {
    emit failed("Unknown error in McWorker::runGreeks");
  }
}

namespace { // helpers locaux

struct McOneResult {
  double mean{0.0};
  double se{0.0};
  std::size_t n{0};
};

McOneResult mc_price_once(const rw::market::MarketData& mkt,
                          const rw::market::Instrument& inst,
                          const rw::config::McConfig& cfg,
                          double sigma,
                          std::atomic<bool>& stopFlag)
{
  using rw::core::RunningStats;
  rw::models::Gbm model({mkt.r, mkt.q, sigma});
  const double T  = inst.T;
  const double df = (T > 0.0 ? std::exp(-mkt.r*T) : 1.0);

  RunningStats acc;
  std::size_t Ncum = 0;
  rw::core::NormalRng rng(cfg.seed);
  const std::size_t nSteps = (cfg.n_steps == 0 ? 1 : cfg.n_steps);
  const double dt = (nSteps <= 1 ? T : (T / static_cast<double>(nSteps)));

  if (T == 0.0) {
    const double pv = (inst.type == rw::market::OptionType::Call)
                      ? rw::payoffs::payoff_call(mkt.S0, inst.K)
                      : rw::payoffs::payoff_put (mkt.S0, inst.K);
    acc.add(pv);
    return {acc.mean(), acc.std_error(), acc.count()};
  }

  while (Ncum < cfg.n_paths_target && !stopFlag.load(std::memory_order_relaxed)) {
    const std::size_t budget = cfg.n_paths_target - Ncum;
    const std::size_t batchN = std::min(cfg.batch_size, budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      double ST = 0.0;
      if (nSteps <= 1) {
        ST = model.sample_ST(mkt.S0, T, rng);
      } else {
        double S = mkt.S0;
        for (std::size_t k=0;k<nSteps;++k) model.step_inplace(S, dt, rng);
        ST = S;
      }
      const double payoff = (inst.type == rw::market::OptionType::Call)
                            ? rw::payoffs::payoff_call(ST, inst.K)
                            : rw::payoffs::payoff_put (ST, inst.K);
      acc.add(df * payoff);
    }

    Ncum += batchN;

    const double se = acc.std_error();
    if (cfg.tolerance > 0.0 && se < cfg.tolerance) break;
  }

  return {acc.mean(), acc.std_error(), acc.count()};
}

} // namespace

void McWorker::runStress(rw::market::MarketData baseMkt,
                         rw::market::Instrument baseInst,
                         rw::config::McConfig cfg,
                         double baseSigma,
                         double dS_rel, double dSigma, double dR, double dT)
{
  stop_.store(false, std::memory_order_relaxed);

  try {
    if (baseMkt.S0 <= 0.0) { emit failed("Invalid baseline S0"); return; }
    if (baseInst.K <= 0.0 || baseInst.T < 0.0) { emit failed("Invalid baseline K/T"); return; }
    if (baseSigma < 0.0) { emit failed("Invalid baseline sigma"); return; }

    qDebug() << "[runStress]"
            << "BASE S0=" << baseMkt.S0 << "r=" << baseMkt.r << "q=" << baseMkt.q
            << "K=" << baseInst.K << "T=" << baseInst.T
            << "sigmaGlobal=" << baseSigma
            << "smileEnabled=" << smileEnabled_ << "vega%=" << vegaStressPct_;
    qDebug() << "[runStress]" << "shocks:"
            << "dS_rel=" << dS_rel << "(=" << dS_rel*100.0 << "%)"
            << "dSigma(abs)=" << dSigma
            << "dR=" << dR
            << "dT=" << dT;
    // 1) Sigma(s) effectives (intègrent Smile + vegaStress)
    const double sigBase = effectiveSigma_(baseMkt, baseInst, baseSigma, /*applyVega=*/false);

    // Construire le scénario stressé (sans muter les structs d’origine)
    const double S0s = baseMkt.S0 * (1.0 + dS_rel);
    const double rs  = baseMkt.r + dR;
    const double qs  = baseMkt.q;                          // pas de choc q ici
    const double Ts  = std::max(0.0, baseInst.T + dT);
    if (S0s <= 0.0) { emit failed("Stressed S0 <= 0"); return; }

    rw::market::MarketData     mktSt(S0s, rs, qs);
    rw::market::Instrument     instSt(baseInst.K, Ts, baseInst.type);

    // sigma stressée = (Smile+vegaStress sur état stressé) + dSigma (absolu)
    double sigSt = effectiveSigma_(mktSt, instSt, baseSigma, /*applyVega=*/true) + dSigma;
    qDebug() << "[runStress]" << "sigmaBase=" << sigBase << "sigmaStress(before +dSigma)="
         << (sigSt - dSigma) << "sigmaStress(final)=" << sigSt;
    if (sigSt <= 0.0) sigSt = 1e-8; 

    // 2) Chrono
    const auto t0 = std::chrono::steady_clock::now();

    rw::core::NormalRng rng(cfg.seed);
    const std::size_t nSteps = (cfg.n_steps == 0 ? 1 : cfg.n_steps);

    double sumBase = 0.0, sumStress = 0.0;
    std::size_t Ncum = 0;

    const double T0 = baseInst.T;
    const double r0 = baseMkt.r, q0 = baseMkt.q;

    // Cas T==0 éventuellement déterministe
    if (T0 == 0.0 && Ts == 0.0) {
      const double pv0 = payoff_from_ST(baseMkt.S0, baseInst);                       // df=1
      const double pvs = payoff_from_ST(S0s,       rw::market::Instrument{baseInst.K, Ts, baseInst.type});
      emit stressFinished(pv0, pvs, (pvs - pv0), 0);
      return;
    }

    while (Ncum < cfg.n_paths_target && !stop_.load(std::memory_order_relaxed)) {
      const std::size_t budget = cfg.n_paths_target - Ncum;
      const std::size_t batchN = std::min(cfg.batch_size, budget);

      for (std::size_t i=0;i<batchN;++i) {
        double ST0=0.0, STs=0.0;

        if (nSteps <= 1) {
          const double z = rng.sample();

          // Exact 1-pas
          const double mu0 = (r0 - q0 - 0.5*sigBase*sigBase)*T0;
          const double mus = (rs - qs - 0.5*sigSt  *sigSt  )*Ts;
          ST0 = baseMkt.S0 * std::exp(mu0 + sigBase*std::sqrt(T0)*z);
          STs = S0s       * std::exp(mus + sigSt  *std::sqrt(Ts)*z);
        } else {
          const double dt0 = (nSteps>1 ? T0/static_cast<double>(nSteps) : T0);
          const double dts = (nSteps>1 ? Ts/static_cast<double>(nSteps) : Ts);

          double S0p = baseMkt.S0;
          double Sps = S0s;

          const double drift0 = (r0 - q0 - 0.5*sigBase*sigBase);
          const double drifts = (rs - qs - 0.5*sigSt  *sigSt);
          for (std::size_t k=0;k<nSteps;++k) {
            const double z = rng.sample();
            if (T0>0.0) S0p *= std::exp(drift0*dt0 + sigBase*std::sqrt(dt0)*z);
            if (Ts>0.0) Sps *= std::exp(drifts*dts + sigSt  *std::sqrt(dts)*z);
          }
          ST0 = S0p; STs = Sps;
        }

        const double df0 = std::exp(-r0*T0);
        const double dfs = std::exp(-rs*Ts);
        const double pv0 = df0 * payoff_from_ST(ST0, baseInst);
        const double pvs = dfs * payoff_from_ST(STs, rw::market::Instrument{baseInst.K, Ts, baseInst.type});

        sumBase   += pv0;
        sumStress += pvs;
      }

      Ncum += batchN;
      // pas d’arrêt par tolérance ici (on compare des moyennes corrélées)
    }

    if (stop_.load(std::memory_order_relaxed)) { emit canceled(); return; }

    const auto t1 = std::chrono::steady_clock::now();
    const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    const double meanBase   = sumBase   / static_cast<double>(Ncum ? Ncum : 1);
    const double meanStress = sumStress / static_cast<double>(Ncum ? Ncum : 1);
    emit stressFinished(meanBase, meanStress, (meanStress - meanBase), ms);

  } catch (const std::exception& e) {
    emit failed(QString::fromUtf8(e.what()));
  } catch (...) {
    emit failed("Unknown error in McWorker::runStress");
  }

}

void McWorker::setSmileMode(bool enabled) {
  smileEnabled_ = enabled;
  qDebug() << "[Worker] setSmileMode =" << enabled;
}

void McWorker::setVegaStressPct(double pct) {
  vegaStressPct_ = pct;
  qDebug() << "[Worker] setVegaStressPct =" << pct;
}

void McWorker::setSmileSurface(std::shared_ptr<const rw::smile::SmileSurface> s) {
  smileSurf_ = std::move(s);
  qDebug() << "[Worker] setSmileSurface ptr?=" << static_cast<bool>(smileSurf_);
}


double McWorker::effectiveSigma_(const rw::market::MarketData& mkt,
                                 const rw::market::Instrument& inst,
                                 double sigmaFallback,
                                 bool applyVega) const {
  double sigma = sigmaFallback;

  double iv = std::numeric_limits<double>::quiet_NaN();
  bool usedSmile = false;
  if (smileEnabled_ && smileSurf_) {
    iv = smileSurf_->iv_at(inst.T, inst.K);
    if (std::isfinite(iv) && iv > 0.0) { sigma = iv; usedSmile = true; }
  }

  bool appliedVega = false;
  if (applyVega && vegaStressPct_ != 0.0) {
    sigma *= (1.0 + vegaStressPct_ / 100.0);
    appliedVega = true;
  }
  const QString ivStr = (std::isfinite(iv) ? QString::number(iv, 'g', 8) : "NA");
  qDebug() << "[Sigma]"
           << "T=" << inst.T << "K=" << inst.K
           << "fallback=" << sigmaFallback
           << "smileEnabled=" << smileEnabled_
           << "iv=" << ivStr
           << "usedSmile=" << usedSmile
           << "applyVega=" << applyVega
           << "vega%=" << vegaStressPct_
           << "appliedVega=" << appliedVega
           << "=> effSigma=" << sigma;
  return sigma;
}


} // namespace gui
