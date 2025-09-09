// src/pricing/greeks.cpp
#include <rw/pricing/greeks.hpp>

#include <rw/payoffs/vanilla.hpp>
#include <rw/core/stats.hpp>

#include <cmath>      // exp, sqrt
#include <algorithm>  // min, max
#include <vector>
#include <chrono>

namespace rw {
namespace pricing {

namespace {

inline double discount(double r, double T) { return std::exp(-r * T); }

inline double payoff_from_ST(double ST, const rw::market::Instrument& inst) {
  using rw::market::OptionType;
  return (inst.type == OptionType::Call)
         ? rw::payoffs::payoff_call(ST, inst.K)
         : rw::payoffs::payoff_put (ST, inst.K);
}

struct Prop {
  double r, q, sigma, T;
  std::size_t nSteps;
  double dt, muT, volT, muDt, volDt;

  Prop(double r_, double q_, double sigma_, double T_, std::size_t steps_)
  : r(r_), q(q_), sigma(std::max(0.0, sigma_)), T(std::max(0.0, T_)), nSteps(steps_) {
    nSteps = (nSteps == 0 ? 1 : nSteps);
    dt = (nSteps <= 1 ? T : (T / static_cast<double>(nSteps)));
    const double s2 = sigma * sigma;
    muT  = (r - q - 0.5 * s2) * T;
    volT = sigma * std::sqrt(T);
    muDt  = (r - q - 0.5 * s2) * dt;
    volDt = sigma * std::sqrt(dt);
  }

  // Évolution avec une séquence de Z (ou un seul si nSteps==1).
  // negate=false => Zk ; negate=true => -Zk.
  double evolve_with(const double S0, const std::vector<double>& Zs, bool negate) const {
    if (T == 0.0) return S0;
    if (nSteps <= 1) {
      double z = Zs[0];
      if (negate) z = -z;
      return S0 * std::exp(muT + volT * z);
    }
    double S = S0;
    for (std::size_t k = 0; k < nSteps; ++k) {
      const double z = negate ? -Zs[k] : Zs[k];
      S *= std::exp(muDt + volDt * z);
    }
    return S;
  }
};

constexpr double Z95 = 1.959963984540054;
inline double half_width_95(double se) { return Z95 * se; }

} // namespace

// ============================================================================
// DELTA (bump relatif S0 : eps)
// ============================================================================
MonteCarloResult GreeksEstimator::delta_brv(const rw::market::MarketData& mkt,
                                            const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params(); // r,q,sigma
  const double r = pm.r, q = pm.q, sigma = pm.sigma;
  const double T = inst.T;

  // T == 0 : delta par FD sur payoff (pas d'actualisation)
  if (T == 0.0) {
    const double eps = gcfg_.bump_rel_S0;
    const double Sp = mkt.S0 * (1.0 + eps);
    const double Sm = mkt.S0 * (1.0 - eps);
    const double Xp = payoff_from_ST(Sp, inst);
    const double Xm = payoff_from_ST(Sm, inst);
    const double delta = (Xp - Xm) / (2.0 * eps * mkt.S0);
    return MonteCarloResult{ delta, 0.0, delta, delta, 0, 0, {} };
  }

  rw::core::NormalRng rng(mccfg_.seed);
  const double eps = gcfg_.bump_rel_S0;
  const double S0p = mkt.S0 * (1.0 + eps);
  const double S0m = mkt.S0 * (1.0 - eps);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);
  const Prop prop(r, q, sigma, T, nSteps);
  const double df = discount(r, T);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;
  const bool useCRN  = gcfg_.use_crn;

  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> ZsA; ZsA.reserve(L);
  std::vector<double> ZsB; ZsB.reserve(L);

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = std::min(mccfg_.batch_size, budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      // 1) Z pour scénario "+", et pour "−" selon CRN
      ZsA.clear();
      for (std::size_t k = 0; k < L; ++k) ZsA.push_back(rng.sample());
      if (useCRN) {
        ZsB = ZsA;
      } else {
        ZsB.clear();
        for (std::size_t k = 0; k < L; ++k) ZsB.push_back(rng.sample());
      }

      // 2) Contribution (avec ou sans antithétique)
      auto contrib_from = [&](const std::vector<double>& Zp,
                              const std::vector<double>& Zm,
                              bool antithetic) -> double {
        // Z : X+ - X- sur S bumpés (même Zp/Zm si CRN)
        const double STp_Z = prop.evolve_with(S0p, Zp, /*neg=*/false);
        const double STm_Z = prop.evolve_with(S0m, Zm, /*neg=*/false);
        const double Xp_Z = df * payoff_from_ST(STp_Z, inst);
        const double Xm_Z = df * payoff_from_ST(STm_Z, inst);
        const double g_Z = (Xp_Z - Xm_Z) / (2.0 * eps * mkt.S0);

        if (!antithetic) return g_Z;

        // -Z : antithétique
        const double STp_nZ = prop.evolve_with(S0p, Zp, /*neg=*/true);
        const double STm_nZ = prop.evolve_with(S0m, Zm, /*neg=*/true);
        const double Xp_nZ = df * payoff_from_ST(STp_nZ, inst);
        const double Xm_nZ = df * payoff_from_ST(STm_nZ, inst);
        const double g_nZ = (Xp_nZ - Xm_nZ) / (2.0 * eps * mkt.S0);

        return 0.5 * (g_Z + g_nZ);
      };

      const double g_i = contrib_from(ZsA, ZsB, useAnti);
      acc.add(g_i);

      // 3) Comptage des chemins physiques consommés
      std::size_t phys = 1;        // ZsA
      if (useAnti)  phys *= 2;     // (Z, -Z)
      if (!useCRN)  phys *= 2;     // ZsB indépendant
      Ncum_phys += phys;
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys;   // chemins physiques consommés // unités statistiques (paires si anti)
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ============================================================================
// VEGA (bump absolu sigma : h) — propagation sous sig+ et sig−
// ============================================================================
MonteCarloResult GreeksEstimator::vega_brv(const rw::market::MarketData& mkt,
                                           const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r = pm.r, q = pm.q, sigma0 = pm.sigma;
  const double T = inst.T;

  if (T == 0.0) {
    return MonteCarloResult{ 0.0, 0.0, 0.0, 0.0, 0, 0, {} };
  }

  const double h = std::max(0.0, gcfg_.bump_abs_sigma);
  const double sig_p = std::max(0.0, sigma0 + h);
  const double sig_m = std::max(0.0, sigma0 - h);

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);
  const Prop prop_p(r, q, sig_p, T, nSteps);
  const Prop prop_m(r, q, sig_m, T, nSteps);
  const double df = discount(r, T);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;
  const bool useCRN  = gcfg_.use_crn;

  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> ZsA; ZsA.reserve(L);
  std::vector<double> ZsB; ZsB.reserve(L);

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = std::min(mccfg_.batch_size, budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      ZsA.clear();
      for (std::size_t k = 0; k < L; ++k) ZsA.push_back(rng.sample());
      if (useCRN) {
        ZsB = ZsA;
      } else {
        ZsB.clear();
        for (std::size_t k = 0; k < L; ++k) ZsB.push_back(rng.sample());
      }

      auto contrib_from = [&](const std::vector<double>& Zp,
                              const std::vector<double>& Zm,
                              bool antithetic) -> double {
        const double STp_Z = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/false);
        const double STm_Z = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/false);
        const double Xp_Z  = df * payoff_from_ST(STp_Z, inst);
        const double Xm_Z  = df * payoff_from_ST(STm_Z, inst);
        const double g_Z   = (Xp_Z - Xm_Z) / (2.0 * h);

        if (!antithetic) return g_Z;

        const double STp_nZ = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/true);
        const double STm_nZ = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/true);
        const double Xp_nZ  = df * payoff_from_ST(STp_nZ, inst);
        const double Xm_nZ  = df * payoff_from_ST(STm_nZ, inst);
        const double g_nZ   = (Xp_nZ - Xm_nZ) / (2.0 * h);

        return 0.5 * (g_Z + g_nZ);
      };

      const double g_i = contrib_from(ZsA, ZsB, useAnti);
      acc.add(g_i);

      std::size_t phys = 1;
      if (useAnti) phys *= 2;
      if (!useCRN) phys *= 2;
      Ncum_phys += phys;
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys;   // chemins physiques consommés
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ============================================================================
// RHO (bump absolu r : dr) — drift (r−q) et discount e^{−rT} changent
// ============================================================================
MonteCarloResult GreeksEstimator::rho_brv(const rw::market::MarketData& mkt,
                                          const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r0 = pm.r, q = pm.q, sigma = pm.sigma;
  const double T  = inst.T;

  const double dr = std::max(0.0, gcfg_.bump_abs_r);

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);

  const double r_p = r0 + dr;
  const double r_m = r0 - dr;

  const Prop prop_p(r_p, q, sigma, T, nSteps);
  const Prop prop_m(r_m, q, sigma, T, nSteps);
  const double df_p = discount(r_p, T);
  const double df_m = discount(r_m, T);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;
  const bool useCRN  = gcfg_.use_crn;

  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> ZsA; ZsA.reserve(L);
  std::vector<double> ZsB; ZsB.reserve(L);

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = std::min(mccfg_.batch_size, budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      ZsA.clear();
      for (std::size_t k = 0; k < L; ++k) ZsA.push_back(rng.sample());
      if (useCRN) { ZsB = ZsA; }
      else { ZsB.clear(); for (std::size_t k = 0; k < L; ++k) ZsB.push_back(rng.sample()); }

      auto contrib_from = [&](const std::vector<double>& Zp,
                              const std::vector<double>& Zm,
                              bool antithetic) -> double {
        // Z
        const double STp_Z = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/false);
        const double STm_Z = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/false);
        const double Xp_Z  = df_p * payoff_from_ST(STp_Z, inst);
        const double Xm_Z  = df_m * payoff_from_ST(STm_Z, inst);
        const double g_Z   = (Xp_Z - Xm_Z) / (2.0 * dr);

        if (!antithetic) return g_Z;

        // -Z
        const double STp_nZ = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/true);
        const double STm_nZ = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/true);
        const double Xp_nZ  = df_p * payoff_from_ST(STp_nZ, inst);
        const double Xm_nZ  = df_m * payoff_from_ST(STm_nZ, inst);
        const double g_nZ   = (Xp_nZ - Xm_nZ) / (2.0 * dr);

        return 0.5 * (g_Z + g_nZ);
      };

      const double g_i = contrib_from(ZsA, ZsB, useAnti);
      acc.add(g_i);

      std::size_t phys = 1;
      if (useAnti) phys *= 2;
      if (!useCRN) phys *= 2;
      Ncum_phys += phys;
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys;   // chemins physiques consommés
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ============================================================================
// THETA (bump absolu T : h) — convention marché: Theta = -dV/dt = -dV/dT
// Central si T>h ; sinon amont : -(V(T+h)-V(T))/h.
// ============================================================================
MonteCarloResult GreeksEstimator::theta_brv(const rw::market::MarketData& mkt,
                                            const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r = pm.r, q = pm.q, sigma = pm.sigma;
  const double T0 = inst.T;

  const double h = std::max(0.0, gcfg_.bump_abs_T);
  const bool central = (T0 > h);

  const double Tp = T0 + h;
  const double Tm = central ? (T0 - h) : T0; // si pas central, Tm=T (diff avant)

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);

  const Prop prop_p(r, q, sigma, Tp, nSteps);
  const Prop prop_m(r, q, sigma, Tm, nSteps);
  const double df_p = discount(r, Tp);
  const double df_m = discount(r, Tm);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;
  const bool useCRN  = gcfg_.use_crn;

  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> ZsA; ZsA.reserve(L);
  std::vector<double> ZsB; ZsB.reserve(L);

  const double denom = central ? (2.0 * h) : h;
  const double sign  = -1.0; // Theta = - dV/dT

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = std::min(mccfg_.batch_size, budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      ZsA.clear();
      for (std::size_t k = 0; k < L; ++k) ZsA.push_back(rng.sample());
      if (useCRN) { ZsB = ZsA; }
      else { ZsB.clear(); for (std::size_t k = 0; k < L; ++k) ZsB.push_back(rng.sample()); }

      auto contrib_from = [&](const std::vector<double>& Zp,
                              const std::vector<double>& Zm,
                              bool antithetic) -> double {
        // Z
        const double STp_Z = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/false);
        const double Xp_Z  = df_p * payoff_from_ST(STp_Z, inst);

        const double STm_Z = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/false);
        const double Xm_Z  = df_m * payoff_from_ST(STm_Z, inst);

        const double g_Z   = sign * (Xp_Z - Xm_Z) / denom;

        if (!antithetic) return g_Z;

        // -Z
        const double STp_nZ = prop_p.evolve_with(mkt.S0, Zp, /*neg=*/true);
        const double Xp_nZ  = df_p * payoff_from_ST(STp_nZ, inst);

        const double STm_nZ = prop_m.evolve_with(mkt.S0, Zm, /*neg=*/true);
        const double Xm_nZ  = df_m * payoff_from_ST(STm_nZ, inst);

        const double g_nZ   = sign * (Xp_nZ - Xm_nZ) / denom;

        return 0.5 * (g_Z + g_nZ);
      };

      const double g_i = contrib_from(ZsA, ZsB, useAnti);
      acc.add(g_i);

      std::size_t phys = 1;
      if (useAnti) phys *= 2;
      if (!useCRN) phys *= 2;
      Ncum_phys += phys;
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys;   // chemins physiques consommés
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ======================== DELTA PATHWISE ========================
MonteCarloResult GreeksEstimator::delta_pathwise(const rw::market::MarketData& mkt,
                                                 const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r = pm.r, q = pm.q, sigma = pm.sigma;
  const double T = inst.T;
  const double S0 = mkt.S0;

  // Cas déterministes
  if (T == 0.0) {
    using rw::market::OptionType;
    double delta0 = 0.0;
    if (inst.type == OptionType::Call) delta0 = (S0 > inst.K) ? 1.0 : 0.0;
    else                               delta0 = (S0 < inst.K) ? -1.0 : 0.0;
    return MonteCarloResult{ delta0, 0.0, delta0, delta0, 0, 0, {} };
  }
  if (sigma == 0.0) {
    const double ST = S0 * std::exp((r - q) * T);
    const double df = std::exp(-r * T);
    using rw::market::OptionType;
    double gi = 0.0;
    if (inst.type == OptionType::Call) gi = (ST > inst.K) ? (df * (ST / S0)) : 0.0;
    else                               gi = (ST < inst.K) ? (-df * (ST / S0)) : 0.0;
    return MonteCarloResult{ gi, 0.0, gi, gi, 0, 0, {} };
  }

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);
  const Prop prop(r, q, sigma, T, nSteps);
  const double df = discount(r, T);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;

  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> Zs; Zs.reserve(L);

  auto pw_sample = [&](const std::vector<double>& Z, bool negate) -> double {
    const double ST = prop.evolve_with(S0, Z, negate);
    using rw::market::OptionType;
    if (inst.type == OptionType::Call) return (ST > inst.K) ? (df * (ST / S0)) : 0.0;
    else                               return (ST < inst.K) ? (-df * (ST / S0)) : 0.0;
  };

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = (mccfg_.batch_size < budget ? mccfg_.batch_size : budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      Zs.clear();
      for (std::size_t k = 0; k < L; ++k) Zs.push_back(rng.sample());

      const double gZ  = pw_sample(Zs, /*neg=*/false);
      double gi = gZ;
      if (useAnti) {
        const double gNZ = pw_sample(Zs, /*neg=*/true);
        gi = 0.5 * (gZ + gNZ);
        Ncum_phys += 2;
      } else {
        Ncum_phys += 1;
      }
      acc.add(gi);
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys; // chemins physiques (comme tes autres funcs)
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ========== GAMMA via central sur Delta-PW (CRN intra-paire) ==========
MonteCarloResult GreeksEstimator::gamma_from_delta_pw(const rw::market::MarketData& mkt,
                                                      const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r = pm.r, q = pm.q, sigma = pm.sigma;
  const double T = inst.T;
  const double S0 = mkt.S0;
  const double eps = (gcfg_.bump_rel_S0 > 0.0 ? gcfg_.bump_rel_S0 : 0.01);
  const double S0p = S0 * (1.0 + eps);
  const double S0m = S0 * (1.0 - eps);

  // Dérivée PW déterministe (T==0 ou sigma==0)
  auto delta_pw_det = [&](double S00)->double {
    if (T == 0.0) {
      using rw::market::OptionType;
      if (inst.type == OptionType::Call) return (S00 > inst.K) ? 1.0 : 0.0;
      else                               return (S00 < inst.K) ? -1.0 : 0.0;
    } else { // sigma==0
      const double ST = S00 * std::exp((r - q) * T);
      const double df = std::exp(-r * T);
      using rw::market::OptionType;
      if (inst.type == OptionType::Call) return (ST > inst.K) ? (df * (ST / S00)) : 0.0;
      else                               return (ST < inst.K) ? (-df * (ST / S00)) : 0.0;
    }
  };
  if (T == 0.0 || sigma == 0.0) {
    const double dP = delta_pw_det(S0p);
    const double dM = delta_pw_det(S0m);
    const double gamma = (dP - dM) / (2.0 * eps * S0);
    return MonteCarloResult{ gamma, 0.0, gamma, gamma, 0, 0, {} };
  }

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);
  const Prop prop(r, q, sigma, T, nSteps);
  const double df = discount(r, T);
  const bool useAnti = gcfg_.use_antithetic;

  auto delta_pw_sample = [&](double S00, const std::vector<double>& Z, bool negate){
    const double ST = prop.evolve_with(S00, Z, negate);
    using rw::market::OptionType;
    if (inst.type == OptionType::Call) return (ST > inst.K) ? (df * (ST / S00)) : 0.0;
    else                               return (ST < inst.K) ? (-df * (ST / S00)) : 0.0;
  };

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const std::size_t L = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> Zs; Zs.reserve(L);

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = (mccfg_.batch_size < budget ? mccfg_.batch_size : budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      Zs.clear();
      for (std::size_t k = 0; k < L; ++k) Zs.push_back(rng.sample());

      // CRN intra-paire: mêmes Z pour S0p et S0m
      const double dP_Z = delta_pw_sample(S0p, Zs, /*neg=*/false);
      const double dM_Z = delta_pw_sample(S0m, Zs, /*neg=*/false);
      double gi = (dP_Z - dM_Z) / (2.0 * eps * S0);

      if (useAnti) {
        const double dP_nZ = delta_pw_sample(S0p, Zs, /*neg=*/true);
        const double dM_nZ = delta_pw_sample(S0m, Zs, /*neg=*/true);
        const double g_nZ  = (dP_nZ - dM_nZ) / (2.0 * eps * S0);
        gi = 0.5 * (gi + g_nZ);
        Ncum_phys += 2;
      } else {
        Ncum_phys += 1;
      }
      acc.add(gi);
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = Ncum_phys;
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}

// ============================================================================
// VEGA LRM (Likelihood Ratio Method)
// ----------------------------------------------------------------------------
// One-step exact:  L = (Z^2 - 1)/sigma
// Multi-step (dt = T/n): L = ( (sum_k sqrt(dt)*Z_k)^2 - T ) / sigma
// Antithétique : L invariant (Z -> -Z), on moyenne les payoffs puis on multiplie par L.
// ============================================================================
MonteCarloResult GreeksEstimator::vega_lrm(const rw::market::MarketData& mkt,
                                           const rw::market::Instrument& inst) const
{
  const auto& pm = model_.params();
  const double r = pm.r, q = pm.q, sigma = pm.sigma;
  const double T = inst.T;
  const double S0 = mkt.S0;

  // Garde-fous : LRM non défini si T<=0 ou sigma<=0
  if (T <= 0.0 || sigma <= 0.0) {
    MonteCarloResult res{};
    res.price = 0.0; res.std_error = 0.0;
    res.ci_low = 0.0; res.ci_high = 0.0;
    res.n_effective = 0; res.elapsed_ms = 0;
    return res;
  }

  rw::core::NormalRng rng(mccfg_.seed);
  const std::size_t nSteps = (gcfg_.n_steps == 0 ? 1 : gcfg_.n_steps);
  const Prop prop(r, q, sigma, T, nSteps);
  const double df = discount(r, T);

  rw::core::RunningStats acc;
  std::vector<ConvergencePoint> log;
  if (log_enabled_) log.reserve((mccfg_.n_paths_target + mccfg_.batch_size - 1) / mccfg_.batch_size);

  const auto t0 = std::chrono::steady_clock::now();

  std::size_t Ncum_phys = 0;
  const bool useAnti = gcfg_.use_antithetic;

  const std::size_t Lsz = (nSteps <= 1 ? 1 : nSteps);
  std::vector<double> Zs; Zs.reserve(Lsz);

  auto payoff_disc = [&](double ST){ return df * payoff_from_ST(ST, inst); };

  while (Ncum_phys < mccfg_.n_paths_target) {
    const std::size_t budget = mccfg_.n_paths_target - Ncum_phys;
    const std::size_t batchN = (mccfg_.batch_size < budget ? mccfg_.batch_size : budget);

    for (std::size_t i = 0; i < batchN; ++i) {
      Zs.clear();

      double contrib = 0.0;
      std::size_t phys = 1;

      if (nSteps <= 1) {
        // ---- 1 PAS ----
        const double z = rng.sample();
        Zs.push_back(z);

        const double STp = prop.evolve_with(S0, Zs, /*neg=*/false);
        const double Xp  = payoff_disc(STp);
        const double Lp  = (z*z - 1.0) / sigma - T * z;      // *** poids correct ***

        if (useAnti) {
          const double STm = prop.evolve_with(S0, Zs, /*neg=*/true);   // -z
          const double Xm  = payoff_disc(STm);
          const double Lm  = (z*z - 1.0) / sigma + T * z;              // *** L(-z) ***
          contrib = 0.5 * (Xp * Lp + Xm * Lm);
          phys = 2;
        } else {
          contrib = Xp * Lp;
          phys = 1;
        }
      } else {
        // ---- MULTI-PAS ----  A = sum sqrt(dt) * Z_k
        double A = 0.0;
        for (std::size_t k = 0; k < Lsz; ++k) {
          const double zk = rng.sample();
          Zs.push_back(zk);
          A += std::sqrt(prop.dt) * zk;
        }

        const double STp = prop.evolve_with(S0, Zs, /*neg=*/false);
        const double Xp  = payoff_disc(STp);
        const double Lp  = ((A*A)/T - 1.0) / sigma - A;               // *** poids correct ***

        if (useAnti) {
          const double STm = prop.evolve_with(S0, Zs, /*neg=*/true);   // Z -> -Z, donc A -> -A
          const double Xm  = payoff_disc(STm);
          const double Lm  = ((A*A)/T - 1.0) / sigma + A;              // *** L(-A) ***
          contrib = 0.5 * (Xp * Lp + Xm * Lm);
          phys = 2;
        } else {
          contrib = Xp * Lp;
          phys = 1;
        }
      }

      acc.add(contrib);
      Ncum_phys += phys;
    }

    if (log_enabled_) log.push_back({ Ncum_phys, acc.mean(), half_width_95(acc.std_error()) });
    if (mccfg_.tolerance > 0.0 && acc.std_error() < mccfg_.tolerance) break;
  }

  const auto t1 = std::chrono::steady_clock::now();

  MonteCarloResult res{};
  res.price        = acc.mean();
  res.std_error    = acc.std_error();
  auto ci          = rw::core::confidence_interval_95(res.price, res.std_error, acc.count());
  res.ci_low       = ci.low;
  res.ci_high      = ci.high;
  res.n_effective  = acc.count();                 // compte statistique (paires si anti)
  res.elapsed_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  if (log_enabled_) res.convergence_log = std::move(log);
  return res;
}



} // namespace pricing
} // namespace rw
