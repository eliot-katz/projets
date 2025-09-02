#pragma once
/**
 * @file gbm.hpp
 * @brief Modèle Black–Scholes / GBM sous mesure risque-neutre Q.
 *
 * # GbmParams
 * - r     : taux sans risque (décimal).
 * - q     : taux de dividende (décimal).
 * - sigma : volatilité (> 0).
 *
 * # Gbm
 * Sous la mesure risque-neutre Q :
 *    dS_t / S_t = (r - q) dt + sigma dW_t
 *
 * Formule exacte pour S_T :
 *    S_T = S0 * exp( (r - q - 0.5*sigma^2) * T + sigma * sqrt(T) * Z )
 * où Z ~ N(0,1).
 *
 * # Méthodes
 * - sample_ST : un pas exact lognormal (spot → maturité).
 * - step_inplace : mise à jour en place d’un spot sur un pas dt.
 *
 * # Tests recommandés
 * - Moments de S_T (espérance et variance) vs formules fermées :
 *     E[S_T] = S0 * exp((r - q)T)
 *     Var[S_T] = E[S_T]^2 * (exp(sigma^2 T) - 1)
 */

#include <cstddef>   // std::size_t
#include "../core/normal_rng.hpp"

namespace rw {
namespace models {

/// @brief Paramètres du modèle GBM (Black–Scholes).
struct GbmParams {
  double r;     ///< Taux sans risque (décimal).
  double q;     ///< Taux de dividende (décimal).
  double sigma; ///< Volatilité (>= 0).
};

/// @brief Modèle Black–Scholes / GBM sous mesure Q.
class Gbm {
public:
  /// @brief Construit un modèle avec paramètres donnés.
  explicit Gbm(GbmParams p);

  /// @return Référence constante aux paramètres du modèle.
  const GbmParams& params() const noexcept { return params_; }

  /**
   * @brief Tire une réalisation de S_T par un pas exact lognormal.
   * @param S0  Spot initial (> 0).
   * @param T   Horizon en années (> 0).
   * @param rng Générateur normal standard (N(0,1)).
   * @return Réalisation de S_T.
   *
   * Formule utilisée :
   * S_T = S0 * exp( (r - q - 0.5*sigma^2)T + sigma*sqrt(T)*Z )
   */
  double sample_ST(double S0, double T, rw::core::NormalRng& rng) const;

  /**
   * @brief Effectue une mise à jour multiplicative en place sur un pas dt.
   * @param S   Spot (modifié en sortie).
   * @param dt  Pas de temps (> 0).
   * @param rng Générateur normal standard (N(0,1)).
   *
   * Formule :
   * S_{t+dt} = S_t * exp( (r - q - 0.5*sigma^2)dt + sigma*sqrt(dt)*Z )
   */
  void step_inplace(double& S, double dt, rw::core::NormalRng& rng) const;

private:
  GbmParams params_;
};

} // namespace models
} // namespace rw
