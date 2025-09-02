#include <rw/pricing/analytic_bs.hpp>

#include <cmath>    // log, exp, sqrt, erf, erfc
#include <limits>   // optionnel: NaN si tu en as besoin

namespace rw {
namespace pricing {

namespace {
// CDF normale standard via erfc pour une meilleure stabilité numérique.
constexpr double INV_SQRT2 = 0.70710678118654752440084436210484903928; // 1/sqrt(2)

inline double Phi(double x) noexcept {
  // Phi(x) = 0.5 * erfc(-x / sqrt(2))
  return 0.5 * std::erfc(-x * INV_SQRT2);
}

inline double pospart(double x) noexcept {
  return (x > 0.0) ? x : 0.0;
}

} // unnamed namespace

double price_call_bs(double S0, double K, double r, double q, double sigma, double T) noexcept {
  // Cas limites
  if (T == 0.0) {
    return pospart(S0 - K); // pas d'actualisation
  }
  if (sigma == 0.0) {
    const double ST = S0 * std::exp((r - q) * T); // trajectoire déterministe
    const double df = std::exp(-r * T);
    return df * pospart(ST - K);
  }
  // Optionnel: cas S0==0 -> call = 0
  if (S0 == 0.0) {
    return 0.0;
  }

  // Calculs généraux, stables
  const double sqrtT     = std::sqrt(T);
  const double sigSqrtT  = sigma * sqrtT;
  const double df        = std::exp(-r * T);
  const double disc_q    = std::exp(-q * T);
  const double logm      = std::log(S0 / K);
  const double muT       = (r - q + 0.5 * sigma * sigma) * T;

  const double d1 = (logm + muT) / sigSqrtT;
  const double d2 = d1 - sigSqrtT;

  // Formule BS
  const double call = disc_q * S0 * Phi(d1) - df * K * Phi(d2);
  return call;
}

double price_put_bs(double S0, double K, double r, double q, double sigma, double T) noexcept {
  // Cas limites
  if (T == 0.0) {
    return pospart(K - S0); // pas d'actualisation
  }
  if (sigma == 0.0) {
    const double ST = S0 * std::exp((r - q) * T); // trajectoire déterministe
    const double df = std::exp(-r * T);
    return df * pospart(K - ST);
  }
  // Optionnel: cas S0==0 -> put ~ K * e^{-rT}
  if (S0 == 0.0) {
    const double df = std::exp(-r * T);
    return df * K;
  }

  // Calculs généraux, stables
  const double sqrtT     = std::sqrt(T);
  const double sigSqrtT  = sigma * sqrtT;
  const double df        = std::exp(-r * T);
  const double disc_q    = std::exp(-q * T);
  const double logm      = std::log(S0 / K);
  const double muT       = (r - q + 0.5 * sigma * sigma) * T;

  const double d1 = (logm + muT) / sigSqrtT;
  const double d2 = d1 - sigSqrtT;

  // Formule BS
  const double put = df * K * Phi(-d2) - disc_q * S0 * Phi(-d1);
  return put;
}

double put_call_parity_gap(double call, double put,
                           double S0, double K, double r, double q, double T) noexcept {
  const double df_r = std::exp(-r * T);
  const double df_q = std::exp(-q * T);
  // gap = call - put - ( S0 e^{-qT} - K e^{-rT} )
  return call - put - (S0 * df_q - K * df_r);
}

} // namespace pricing
} // namespace rw
