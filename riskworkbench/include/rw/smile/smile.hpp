#pragma once
#include <vector>
#include <cstddef>

namespace rw::smile {

/**
 * Smile 1D à maturité fixe T.
 * x = ln(K/S0), y = IV (σ) en ce x.
 * Interpolation: PCHIP (monotone cubic Hermite). Extrapolation: linéaire via pente d’extrémité.
 */
struct Smile1D {
  double T{0.0};           // maturité de la slice
  double S0{1.0};          // spot de référence pour x=ln(K/S0)

  // Données (croissantes en x)
  std::vector<double> x;   // x_i = ln(K_i / S0)
  std::vector<double> y;   // y_i = iv(x_i)

  // Tangentes (slopes) PCHIP aux noeuds, calculées par build()
  std::vector<double> d;   // même taille que x,y

  // Calcule d[] (slopes) selon l’algo PCHIP (Fritsch–Carlson).
  // Hypothèses: x strictement croissante, |x|>=2 et |x|==|y|.
  void build();

  // Évalue l’IV au strike K (via x=ln(K/S0)). Renvoie NaN si données insuffisantes.
  double iv_at(double K) const;

  // Helper: évalue à un x donné (si tu veux déjà être en log-moneyness).
  double iv_at_x(double xq) const;
};

/**
 * Surface de smiles: une slice PCHIP par maturité.
 * Interpolation en temps: linéaire sur la variance implicite w = σ² T.
 */
struct SmileSurface {
  double S0{1.0}, r{0.0}, q{0.0};
  std::vector<Smile1D> slices; // trié par T croissante

  // σ(T,K) = sqrt( w(T,K) / T ), où w(T,K) = interp linéaire entre slices en w.
  double iv_at(double T, double K) const;
};

} // namespace rw::smile
