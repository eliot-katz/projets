#pragma once
/**
 * @file instrument.hpp
 * @brief Description d’une option vanille européenne.
 *
 * # Contenu
 * - Type d’option : Call ou Put.
 * - Strike K (>0).
 * - Maturité T (>0, en années fractionnelles).
 *
 * # Domaine valide
 * - K > 0
 * - T > 0
 *
 * # Convention
 * - Temps T en années fractionnelles (ex : 0.5 = 6 mois).
 * - Cohérence avec la mesure risque-neutre utilisée en Black–Scholes.
 */

#include <stdexcept> // std::invalid_argument

namespace rw {
namespace market {

/// @brief Type d’option vanille (Call ou Put).
enum class OptionType {
  Call, ///< Droit d’acheter l’actif sous-jacent.
  Put   ///< Droit de vendre l’actif sous-jacent.
};

/// @brief Instrument de marché : option vanille européenne.
/// @details Strike et maturité doivent être strictement positifs.
struct Instrument {
public:
  const double K;        ///< Strike (> 0).
  const double T;        ///< Maturité en années (> 0).
  const OptionType type; ///< Call ou Put.

  /// @brief Construit un instrument valide.
  /// @param K    Strike (doit être > 0).
  /// @param T    Maturité en années (doit être > 0).
  /// @param type Type d’option (Call ou Put).
  /// @throws std::invalid_argument si K <= 0 ou T <= 0.
  Instrument(double K, double T, OptionType type)
      : K(K), T(T), type(type) {
    if (K <= 0.0) {
      throw std::invalid_argument("Instrument: K must be > 0");
    }
    if (T <= 0.0) {
      throw std::invalid_argument("Instrument: T must be > 0");
    }
  }
};

} // namespace market
} // namespace rw
