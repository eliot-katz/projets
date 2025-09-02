#pragma once
/**
 * @file vanilla.hpp
 * @brief Payoffs d’options vanille européennes.
 *
 * # Définitions
 * - Call européen : max(S_T - K, 0).
 * - Put européen  : max(K - S_T, 0).
 *
 * # Propriétés analytiques
 * - Continuité : les payoffs sont continus en S_T.
 * - Régularité : non dérivables en S_T = K (point de non-lissité).
 *
 * Ces propriétés sont importantes pour le choix d’estimateurs des Greeks :
 * - Delta d’un call/put est discontinu en K.
 * - Les méthodes de différentiation de Monte Carlo doivent en tenir compte.
 */

#include <algorithm> // std::max

namespace rw {
namespace payoffs {

/// @brief Payoff d’un call européen.
/// @param ST Spot à maturité.
/// @param K  Strike (>0).
/// @return max(ST - K, 0).
inline double payoff_call(double ST, double K) noexcept {
  return std::max(ST - K, 0.0);
}

/// @brief Payoff d’un put européen.
/// @param ST Spot à maturité.
/// @param K  Strike (>0).
/// @return max(K - ST, 0).
inline double payoff_put(double ST, double K) noexcept {
  return std::max(K - ST, 0.0);
}

} // namespace payoffs
} // namespace rw
