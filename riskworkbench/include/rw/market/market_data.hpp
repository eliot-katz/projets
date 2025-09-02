#pragma once
/**
 * @file market_data.hpp
 * @brief Données de marché minimales pour Black–Scholes à paramètres constants.
 *
 * # Contenu
 * - Spot initial S0 (> 0).
 * - Taux sans risque r (décimal, ex : 0.05 pour 5 %).
 * - Taux de dividende q (décimal).
 *
 * # Domaines valides
 * - S0 > 0 obligatoire.
 * - Les taux peuvent être négatifs (cas de marchés réels).
 * - La volatilité sigma >= 0 sera ajoutée plus tard.
 *
 * # Unités
 * - Temps T exprimé en années fractionnelles (ex : 0.5 = 6 mois).
 * - Taux continus annualisés en décimal.
 *
 * # Extensions prévues
 * - Hooks pour courbes de taux et dividendes.
 */

#include <stdexcept> // std::invalid_argument

namespace rw {
namespace market {

/**
 * @brief Données de marché minimales (BS constant).
 *
 * Immuables après construction.
 */
struct MarketData {
public:
  const double S0; ///< Spot initial (> 0).
  const double r;  ///< Taux sans risque (décimal).
  const double q;  ///< Taux de dividende (décimal).

  /// @brief Construit une structure de données de marché valide.
  /// @param S0 Spot initial (doit être > 0).
  /// @param r  Taux sans risque (décimal, peut être négatif).
  /// @param q  Taux de dividende (décimal).
  /// @throws std::invalid_argument si S0 <= 0.
  MarketData(double S0, double r, double q) : S0(S0), r(r), q(q) {
    if (S0 <= 0.0) {
      throw std::invalid_argument("MarketData: S0 must be > 0");
    }
  }
};

} // namespace market
} // namespace rw
