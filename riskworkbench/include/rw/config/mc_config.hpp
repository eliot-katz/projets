#pragma once
/**
 * @file mc_config.hpp
 * @brief Configuration standard d’un run Monte Carlo.
 *
 * # Contenu
 * - n_paths_target : nombre maximal de trajectoires simulées **(chemins physiques)**.
 * - batch_size     : nombre de trajectoires simulées par lot (≥ 1).
 * - tolerance      : précision cible (erreur standard). Si <= 0 ⇒ désactivé.
 * - n_steps        : pas de discrétisation (1 = S_T direct ; >1 ⇒ calendrier uniforme).
 * - seed           : graine du RNG.
 * - use_antithetic : si true, utilise des **paires antithétiques** (Z, -Z).
 *
 * # Control Variates (CV)
 * - use_control_variate : active le schéma de variables de contrôle.
 * - cv_mode             : SameInstrument (CV = même instrument) ou CustomInstrument.
 * - cv_instrument       : instrument de contrôle si CustomInstrument (T typiquement = T cible).
 * - cv_beta_mode        : Pilot (phase pilote pour estimer β) ou Fixed (β fixé).
 * - cv_beta_fixed       : valeur de β si Fixed.
 * - cv_pilot_paths      : budget de **chemins physiques** pour la phase pilote.
 *
 * # Logique d’arrêt
 * - Si tolerance > 0 : arrêt anticipé lorsque l’erreur standard < tolerance,
 *   ou après avoir simulé n_paths_target **chemins physiques** (pilot + phase principale).
 * - Si tolerance <= 0 : toujours simuler n_paths_target **chemins physiques** (pilot inclus).
 */

#include <cstddef>   // std::size_t
#include <cstdint>   // std::uint64_t
#include <rw/market/instrument.hpp> // pour cv_instrument

namespace rw {
namespace config {

enum class CvMode { SameInstrument, CustomInstrument };
enum class BetaMode { Pilot, Fixed };

/// @brief Configuration d’un run Monte Carlo.
struct McConfig {
  std::size_t  n_paths_target;   ///< Nombre maximal de trajectoires (chemins physiques).
  std::size_t  batch_size;       ///< Taille d’un lot (en chemins physiques).
  double       tolerance;        ///< Précision cible (<=0 désactive).
  std::size_t  n_steps;          ///< Nombre de pas (1 = S_T direct).
  std::uint64_t seed;            ///< Graine du RNG.

  bool use_antithetic;           ///< Active les paires antithétiques (Z, -Z).

  // --- Control Variates ---
  bool use_control_variate;            ///< Active la variable de contrôle.
  CvMode cv_mode;                      ///< SameInstrument ou CustomInstrument.
  rw::market::Instrument cv_instrument;///< Instrument de contrôle (si CustomInstrument).
  BetaMode cv_beta_mode;               ///< Pilot (par défaut) ou Fixed.
  double cv_beta_fixed;                ///< Valeur de β si mode Fixed.
  std::size_t cv_pilot_paths;          ///< Budget pilote (chemins physiques).

  /// @brief Construit une configuration avec valeurs par défaut.
  McConfig(std::size_t n_paths_target = 1'000'000,
           std::size_t batch_size = 100'000,
           double      tolerance = -1.0,
           std::size_t n_steps = 1,
           std::uint64_t seed = 42ULL,
           bool use_antithetic = false,
           bool use_control_variate = false,
           CvMode cv_mode = CvMode::SameInstrument,
           rw::market::Instrument cv_instrument = rw::market::Instrument{1.0, 1.0, rw::market::OptionType::Call},
           BetaMode cv_beta_mode = BetaMode::Pilot,
           double cv_beta_fixed = 0.0,
           std::size_t cv_pilot_paths = 20000) noexcept
      : n_paths_target(n_paths_target),
        batch_size(batch_size),
        tolerance(tolerance),
        n_steps(n_steps),
        seed(seed),
        use_antithetic(use_antithetic),
        use_control_variate(use_control_variate),
        cv_mode(cv_mode),
        cv_instrument(cv_instrument),
        cv_beta_mode(cv_beta_mode),
        cv_beta_fixed(cv_beta_fixed),
        cv_pilot_paths(cv_pilot_paths) {}
};

} // namespace config
} // namespace rw
