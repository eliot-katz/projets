#pragma once
/**
 * @file mc_config.hpp
 * @brief Configuration standard d’un run Monte Carlo.
 *
 * # Contenu
 * - n_paths_target : nombre maximal de trajectoires simulées.
 * - batch_size     : nombre de trajectoires simulées par lot (≥ 1).
 * - tolerance      : précision cible (erreur standard). Si <= 0 ⇒ désactivé.
 * - n_steps        : pas de discrétisation (1 = S_T direct ; >1 ⇒ calendrier uniforme).
 * - seed           : graine du RNG.
 *
 * # Valeurs par défaut
 * - n_paths_target = 1'000'000
 * - batch_size     = 100'000
 * - tolerance      = -1.0 (désactivé par défaut)
 * - n_steps        = 1
 * - seed           = 42
 *
 * # Recommandations
 * - batch_size ~ 10^4 à 10^5 : compromis entre performance et précision statistique.
 * - Un batch_size trop petit accroît la variance de l’estimateur intermédiaire.
 *
 * # Logique d’arrêt
 * - Si tolerance > 0 : arrêt anticipé lorsque l’erreur standard < tolerance,
 *   ou après avoir simulé n_paths_target.
 * - Si tolerance <= 0 : toujours simuler n_paths_target trajectoires.
 */

#include <cstddef>   // std::size_t
#include <cstdint>   // std::uint64_t

namespace rw {
namespace config {

/// @brief Configuration d’un run Monte Carlo.
struct McConfig {
  std::size_t n_paths_target; ///< Nombre maximal de trajectoires.
  std::size_t batch_size;     ///< Taille d’un lot de simulation.
  double tolerance;           ///< Précision cible (<=0 désactive).
  std::size_t n_steps;        ///< Nombre de pas (1 = S_T direct).
  std::uint64_t seed;         ///< Graine du RNG.

  /// @brief Construit une configuration avec valeurs par défaut.
  McConfig(std::size_t n_paths_target = 1'000'000,
           std::size_t batch_size = 100'000,
           double tolerance = -1.0,
           std::size_t n_steps = 1,
           std::uint64_t seed = 42ULL) noexcept
      : n_paths_target(n_paths_target),
        batch_size(batch_size),
        tolerance(tolerance),
        n_steps(n_steps),
        seed(seed) {}
};

} // namespace config
} // namespace rw
