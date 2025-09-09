#pragma once
/**
 * @file greeks_config.hpp
 * @brief Configuration des estimations de Greeks par Monte Carlo (bump & reprice).
 *
 * # Objet
 * Décrit le Greek ciblé et les tailles de perturbations ("bumps") à appliquer
 * lors d'un schéma de type bump-and-reprice. La logique d'échantillonnage
 * (CRN, antithétiques, n_steps) est alignée sur le reste de la lib.
 *
 * # Conventions
 * - T en années fractionnelles (ex : 1/365 ≈ 1 jour).
 * - Les "bumps" suivants sont utilisés par défaut :
 *     * S0 : bump relatif (par ex. 1% pour Delta)
 *     * sigma : bump absolu (par ex. +0.01 = +1 vol point)
 *     * r : bump absolu (par ex. +0.0001 = +1 bp)
 *     * T : bump absolu (par ex. +1/365 ≈ +1 jour ; Theta ≈ -dV/dt)
 *
 * # Variance reduction
 * - use_crn : si true, utilise les mêmes tirages (Common Random Numbers) entre
 *   le prix "base" et les prix "bumpés" pour réduire la variance de l'estimateur.
 * - use_antithetic : si true, réutilise la logique antithétique (Z, -Z).
 *
 * # Remarques
 * - Cette config ne duplique pas n_paths_target / batch / tolerance / seed.
 *   On réutilise votre rw::config::McConfig pour ces paramètres d'exécution.
 */

#include <cstddef> // std::size_t

namespace rw {
namespace config {

/// @brief Greeks pris en charge par les schémas bump-and-reprice.
enum class Greek {
  Delta,  ///< dV/dS0
  Vega,   ///< dV/dsigma
  Gamma,  ///< d^2V/dS0^2 (schéma central recommandé)
  Rho,    ///< dV/dr
  Theta   ///< -dV/dt (attention : bump sur T, en années)
};

/// @brief Configuration d’un calcul de Greek (bump & reprice).
struct GreeksConfig {
  Greek greek;                ///< Greek visé.

  // Tailles de bump par défaut
  double bump_rel_S0    = 0.01;      ///< Bump relatif sur S0 (ex: 1% pour Delta/Gamma).
  double bump_abs_sigma = 0.01;      ///< Bump absolu sur sigma (ex: +1 vol point pour Vega).
  double bump_abs_r     = 0.0001;    ///< Bump absolu sur r (ex: +1 bp pour Rho).
  double bump_abs_T     = 1.0 / 365; ///< Bump absolu sur T (ex: +1 jour pour Theta).

  // Variance reduction / discrétisation
  bool use_crn         = true;   ///< Common Random Numbers entre base et bumps.
  bool use_antithetic  = false;  ///< Réutilise Z/-Z pour chaque run (base et bumps).
  std::size_t n_steps  = 1;      ///< 1 = pas exact vers T ; >1 = calendrier uniforme.

  // --- Pathwise Delta ---
  bool use_pathwise     = false;  // active collecte PW (Delta). Gamma utilise central sur Delta-PW

  /// @brief Constructeur explicite avec choix du Greek (le reste a des défauts sensés).
  explicit GreeksConfig(Greek g) : greek(g) {}

  /// @brief Constructeur complet si besoin de tout préciser.
  GreeksConfig(Greek g,
               double bumpRelS0,
               double bumpAbsSigma,
               double bumpAbsR,
               double bumpAbsT,
               bool   useCrn,
               bool   useAnti,
               std::size_t nSteps)
  : greek(g),
    bump_rel_S0(bumpRelS0),
    bump_abs_sigma(bumpAbsSigma),
    bump_abs_r(bumpAbsR),
    bump_abs_T(bumpAbsT),
    use_crn(useCrn),
    use_antithetic(useAnti),
    n_steps(nSteps) {}
};

} // namespace config
} // namespace rw
