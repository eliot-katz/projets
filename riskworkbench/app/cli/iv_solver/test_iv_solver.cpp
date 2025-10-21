#include "rw/pricing/implied_vol.hpp"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

// BS de référence (pour générer les prix)
static double norm_cdf(double x){ return 0.5 * std::erfc(-x/std::sqrt(2.0)); }
static double bs_call(double S, double K, double r, double q, double T, double sigma) {
  const double df_r = std::exp(-r*T), df_q = std::exp(-q*T);
  if (T<=0.0 || sigma<=0.0) return std::max(0.0, S*df_q - K*df_r);
  const double vol = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T) / vol;
  const double d2 = d1 - vol;
  return S*df_q*norm_cdf(d1) - K*df_r*norm_cdf(d2);
}
static double bs_put(double S, double K, double r, double q, double T, double sigma) {
  const double df_r = std::exp(-r*T), df_q = std::exp(-q*T);
  if (T<=0.0 || sigma<=0.0) return std::max(0.0, K*df_r - S*df_q);
  const double vol = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T) / vol;
  const double d2 = d1 - vol;
  return K*df_r*norm_cdf(-d2) - S*df_q*norm_cdf(-d1);
}

int main() {
  const double S0 = 100.0, r = 0.01, q = 0.00, T = 0.75;
  const double SIG_TRUE = 0.30;

  // Grille de strikes (moneyness varié)
  std::vector<double> Ks {60, 70, 80, 90, 100, 110, 120, 140, 160};

  // 1) Exactitude : prix synthétiques → retrouve σ à 1e-6 près (calls & puts)
  for (double K : Ks) {
    const double C = bs_call(S0, K, r, q, T, SIG_TRUE);
    const auto ivC = rw::bs::implied_vol_cp(S0, K, r, q, T, C, /*is_call=*/true);
    assert(ivC.converged);
    assert(std::abs(ivC.sigma - SIG_TRUE) < 1e-6);

    const double P = bs_put(S0, K, r, q, T, SIG_TRUE);
    const auto ivP = rw::bs::implied_vol_cp(S0, K, r, q, T, P, /*is_call=*/false);
    assert(ivP.converged);
    assert(std::abs(ivP.sigma - SIG_TRUE) < 1e-6);
  }

  // 2) Contrôle “pas de zig-zag” : IV ~ constante sur K si vol vraie constante
  std::vector<double> ivs;
  ivs.reserve(Ks.size());
  for (double K : Ks) {
    const double C = bs_call(S0, K, r, q, T, SIG_TRUE);
    ivs.push_back(rw::bs::implied_vol(S0, K, r, q, T, C).sigma);
  }
  // max écart à SIG_TRUE
  double devmax = 0.0;
  for (double v : ivs) devmax = std::max(devmax, std::abs(v - SIG_TRUE));
  // tolérance un peu plus large à cause de l’échelle prix/vega (mais très serrée)
  assert(devmax < 5e-7);

  std::cout << "IV solver OK. Max dev=" << devmax << "\n";
  return 0;
}
