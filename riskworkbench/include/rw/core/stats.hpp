#pragma once
/**
 * @file stats.hpp
 * @brief Accumulateurs statistiques en streaming (Welford) + IC 95 %.
 *
 * - Algorithme de Welford : stable numériquement, une passe.
 * - Variance : échantillon (diviseur n-1).
 * - Erreur standard : sqrt(variance / n).
 * - IC 95 % : mean ± z * std_error, avec z ≈ 1.9599639845 (normale).
 *
 * Remarques :
 * - Pas de dépendance lourde dans ce header (uniquement <cstddef>).
 * - Les fonctions nécessitant sqrt sont déclarées ici et peuvent
 *   être définies dans un .cpp qui inclut <cmath>.
 * - Comportement aux petits n :
 *   - n == 0 : mean()=NaN? (non défini ici), variance()=NaN?, std_error()=NaN?
 *   - n == 1 : variance()=0, std_error()=0.
 *   Tu peux choisir la politique dans l’implémentation (NaN vs 0), mais
 *   documente-la dans le .cpp pour tes tests.
 */

#include <cstddef> // std::size_t

namespace rw {
namespace core {

struct RunningStats {
public:
  /// @brief Initialise les accumulateurs (n=0, mean=0, M2=0).
  RunningStats() noexcept;

  /// @brief Ajoute un échantillon.
  void add(double x) noexcept;

  /// @return Nombre d’échantillons vus.
  std::size_t count() const noexcept;

  /// @return Moyenne courante.
  double mean() const noexcept;

  /// @return Variance d'échantillon (diviseur n-1).
  /// @note Si n < 2, retourne NaN (voir implémentation).
  [[nodiscard]] double variance() const noexcept;

  /// @return Erreur standard de la moyenne : sqrt(variance / n).
  /// @note Si n == 0, retourne NaN (voir implémentation).
  [[nodiscard]] double std_error() const noexcept;

private:
  std::size_t n_{0};
  double mean_{0.0};
  double m2_{0.0}; // somme des carrés des écarts à la moyenne (Welford)
};

/// @brief Intervalle de confiance 95 % pour la moyenne (approx. normale).
struct ConfidenceInterval {
  double low;
  double high;
};

/// @brief Calcule mean ± z * std_error (z ≈ 1.9599639845).
/// @param mean       moyenne estimée
/// @param std_error  erreur standard de la moyenne
/// @param n          effectif (utilisé si tu veux adapter z pour petits n via t de Student)
/// @return [low, high]
[[nodiscard]] ConfidenceInterval confidence_interval_95(double mean,
                                          double std_error,
                                          std::size_t n) noexcept;

} // namespace core
} // namespace rw
