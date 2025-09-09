#pragma once
/**
 * @file greeks.hpp
 * @brief Estimation de Greeks par Monte Carlo (bump & reprice en un seul passage, CRN/antithétiques).
 *
 * # Principe
 * - Un **seul passage** sur les tirages : pour chaque trajectoire simulée, on calcule
 *   les deux payoffs "bumpés" (±) **avec les mêmes Gaussiennes** (CRN), puis on
 *   accumule directement l'estimateur **différentiel** (ex. différence centrale).
 * - Si activé, on applique aussi la variance réduction **antithétique** (Z, -Z)
 *   de manière cohérente aux deux payoffs bumpés.
 *
 * # Estimateurs implémentés
 * - Delta (central, bump relatif S0) :
 *     Delta ≈ [ V(S0*(1+h)) − V(S0*(1−h)) ] / (2 h S0)
 * - Vega  (central, bump absolu sigma) :
 *     Vega  ≈ [ V(sigma+h)  − V(sigma−h)  ] / (2 h)
 *
 * # Sortie
 * - Retourne un `MonteCarloResult` dont les champs représentent l’estimation du **Greek**
 *   (et non un prix) : `price` = estimateur, `std_error` = erreur standard de l’estimateur,
 *   `ci_low/ci_high` = IC 95 % sur le Greek, `n_effective` = nb d’unités statistiques
 *   (paires si antithétique), `elapsed_ms` = durée, `convergence_log` optionnel.
 *
 * # Conventions
 * - T en années fractionnelles. r, q en décimal.
 * - Les tailles de bump et les options CRN/antithétiques/n_steps viennent de `GreeksConfig`.
 * - Les paramètres d’exécution (n_paths_target, batch, tolerance, seed) viennent de `McConfig`.
 *
 * # Remarques d’implémentation (dans le .cpp)
 * - Aucun second pricer n’est lancé : on fait **un seul passage**, en calculant
 *   pour chaque chemin les deux valuations bumpées et leur contribution différentielle.
 * - En mode antithétique, l’**unité statistique** est la **paire** : on moyenne
 *   la contribution + et −, puis on accumule.
 * - On réutilise `models::Gbm` pour sampler (exact ou multipas via `n_steps`).
 */

#include <cstddef> // std::size_t

#include <rw/models/gbm.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/config/greeks_config.hpp>
#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/pricing/mc_pricer.hpp> // pour MonteCarloResult & ConvergencePoint

namespace rw {
namespace pricing {

/// @brief Estimateur de Greeks par MC (bump & reprice, CRN/antithétiques, un seul passage).
class GreeksEstimator {
public:
  /**
   * @brief Construit un estimateur de Greeks.
   * @param model  Modèle GBM (r, q, sigma) sous Q.
   * @param gcfg   Configuration du Greek (type, bumps, CRN/anti, n_steps).
   * @param mccfg  Config Monte Carlo (n_paths_target, batch, tolerance, seed).
   */
  explicit GreeksEstimator(rw::models::Gbm model,
                           rw::config::GreeksConfig gcfg,
                           rw::config::McConfig mccfg) noexcept
    : model_(model), gcfg_(gcfg), mccfg_(mccfg), log_enabled_(false) {}

  /// @brief Active/désactive le journal de convergence (par lot).
  void enable_convergence_log(bool on) noexcept { log_enabled_ = on; }

  /**
   * @brief Delta par bump-and-reprice central (CRN, anti optionnel).
   * @details
   *   h = gcfg_.bump_rel_S0 ;
   *   Delta ≈ [ V(S0*(1+h)) − V(S0*(1−h)) ] / (2 h S0).
   *   Un seul passage : pour chaque chemin, on calcule les deux payoffs bumpés
   *   avec les mêmes tirages (et leurs antithétiques si activés), puis on
   *   accumule directement la **contribution différentielle**.
   * @return MonteCarloResult (estimateur du Greek + IC 95 % + meta).
   */
  [[nodiscard]] MonteCarloResult
  delta_brv(const rw::market::MarketData& mkt,
            const rw::market::Instrument& inst) const;

  /**
   * @brief Vega par bump-and-reprice central (CRN, anti optionnel).
   * @details
   *   h = gcfg_.bump_abs_sigma ;
   *   Vega ≈ [ V(sigma+h) − V(sigma−h) ] / (2 h).
   *   Même logique : un seul passage, contributions différentielles par chemin.
   *   NB : clamp interne si sigma±h < 0 (éviter sigma < 0).
   * @return MonteCarloResult (estimateur du Greek + IC 95 % + meta).
   */
  [[nodiscard]] MonteCarloResult
  vega_brv(const rw::market::MarketData& mkt,
           const rw::market::Instrument& inst) const;
  
  [[nodiscard]] MonteCarloResult
  rho_brv(const rw::market::MarketData& mkt,
          const rw::market::Instrument& inst) const;

  [[nodiscard]] MonteCarloResult
  theta_brv(const rw::market::MarketData& mkt,
            const rw::market::Instrument& inst) const;

  // Delta Pathwise (PW)
  [[nodiscard]] MonteCarloResult 
  delta_pathwise(const rw::market::MarketData& mkt,
                 const rw::market::Instrument& inst) const;

  // Gamma via central sur Delta-PW, CRN (mêmes Z pour S0+ et S0-)
  [[nodiscard]] MonteCarloResult 
  gamma_from_delta_pw(const rw::market::MarketData& mkt,
                      const rw::market::Instrument& inst) const;

  // Vega par Likelihood Ratio Method (LRM)
  [[nodiscard]] MonteCarloResult 
  vega_lrm(const rw::market::MarketData& mkt,
                      const rw::market::Instrument& inst) const;

  // (Plus tard)
  // [[nodiscard]] MonteCarloResult gamma_brv (const market::MarketData&, const market::Instrument&) const;
  
  /// @return Accès lecture seule.
  const rw::models::Gbm&            model()  const noexcept { return model_;  }
  const rw::config::GreeksConfig&   gcfg()   const noexcept { return gcfg_;   }
  const rw::config::McConfig&       mccfg()  const noexcept { return mccfg_;  }

private:
  rw::models::Gbm          model_;
  rw::config::GreeksConfig gcfg_;
  rw::config::McConfig     mccfg_;
  bool                     log_enabled_;
};

} // namespace pricing
} // namespace rw
