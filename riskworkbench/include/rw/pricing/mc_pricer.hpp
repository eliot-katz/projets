#pragma once
/**
 * @file mc_pricer.hpp
 * @brief Orchestrateur de pricing Monte Carlo générique pour un européen (BS phase 1).
 *
 * # Principe
 * - Simulation sous la mesure risque-neutre Q avec un modèle GBM (Black–Scholes).
 * - Actualisation par exp(-r * T).
 * - Discrétisation : `n_steps` pas uniformes si `n_steps > 1`, sinon pas exact direct vers S_T.
 * - Pas de stockage de chemins (streaming) : accumulation en ligne via RunningStats (variance empirique).
 *
 * # Convergence & arrêt
 * - Deux modes :
 *   1) `tolerance > 0` : arrêt anticipé dès que l’erreur standard < tolerance,
 *      ou au plus tard à `n_paths_target`.
 *   2) `tolerance <= 0` : toujours simuler `n_paths_target`.
 * - Retourne toujours un IC 95 % cohérent : moyenne ± 1.9599 * std_error.
 *
 * # Journal de convergence (optionnel)
 * - Si activé via `enable_convergence_log(true)`, enregistre des points
 *   (n_cum, estimate, half_width_95) à cadence lot/batch.
 *
 * # Unités et conventions
 * - T en années fractionnelles, r et q en décimal (peuvent être négatifs).
 */

#include <vector>
#include <cstddef>

// Headers "maison"
/*
#include "../models/gbm.hpp"
#include "../payoffs/vanilla.hpp"
#include "../core/stats.hpp"
#include "../config/mc_config.hpp"
#include "../market/market_data.hpp"
#include "../market/instrument.hpp"
*/

#include <rw/models/gbm.hpp>
#include <rw/payoffs/vanilla.hpp>
#include <rw/core/stats.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>


namespace rw {
namespace pricing {

/// @brief Un point du journal de convergence (après un cumul de trajectoires).
struct ConvergencePoint {
  std::size_t n_cum;     ///< Nombre cumulé de trajectoires simulées.
  double      estimate;  ///< Estimation courante du prix (moyenne actualisée).
  double      half_width_95; ///< Demi-largeur de l’IC 95 % (≈ 1.9599 * std_error).
};

/// @brief Résultat d’un run Monte Carlo.
struct MonteCarloResult {
  double price;            ///< Estimation finale du prix.
  double std_error;        ///< Erreur standard empirique.
  double ci_low;           ///< Borne basse de l’IC 95 %.
  double ci_high;          ///< Borne haute de l’IC 95 %.
  std::size_t n_effective; ///< Nombre effectif de trajectoires simulées.
  long long elapsed_ms;    ///< Durée de la simulation (millisecondes).
  std::vector<ConvergencePoint> convergence_log; ///< Journal (peut être vide si désactivé).
};

/**
 * @brief Monte Carlo pricer (européen) orchestré sous GBM.
 *
 * Implémentation prévue :
 * - Tirage exact lognormal si `n_steps == 1` via `Gbm::sample_ST`.
 * - Sinon, propagation multipas uniforme via `Gbm::step_inplace`.
 * - Payoffs européens via `rw::payoffs::payoff_call/put` selon `Instrument::type`.
 * - Actualisation par `std::exp(-r * T)` ; pas de stockage des chemins.
 */
class McPricer {
public:
  /// @brief Construit un pricer avec un modèle et une configuration.
  explicit McPricer(rw::models::Gbm model, rw::config::McConfig cfg) noexcept
    : model_(model), cfg_(cfg), log_enabled_(false) {}

  /// @brief Active/désactive le journal de convergence (par lot).
  void enable_convergence_log(bool on) noexcept { log_enabled_ = on; }

  /**
   * @brief Prix d’une option européenne (call/put) sous GBM.
   * @param mkt  Données de marché (S0, r, q).
   * @param inst Instrument (K, T, type).
   * @return MonteCarloResult avec IC 95 % cohérent et méta-infos de run.
   *
   * Détails d’implémentation attendus dans le .cpp :
   * - Utiliser `cfg_.n_steps` pour le calendrier :
   *     * 1       ⇒ pas exact direct vers S_T (lognormal).
   *     * >1      ⇒ calendrier uniforme de pas `dt = T / n_steps`.
   * - Accumuler en streaming via `core::RunningStats` le payoff actualisé.
   * - Construire l’IC 95 % avec `core::confidence_interval_95`.
   * - Respecter la logique d’arrêt : `tolerance` vs `n_paths_target`.
   */
  MonteCarloResult price_european(const rw::market::MarketData& mkt,
                                  const rw::market::Instrument& inst) const;

  /// @return Accès en lecture aux paramètres du modèle (copie profonde non nécessaire ici).
  const rw::models::Gbm& model() const noexcept { return model_; }

  /// @return Configuration Monte Carlo utilisée.
  const rw::config::McConfig& config() const noexcept { return cfg_; }

  // (Optionnel, plus tard) : support d’un calendrier non uniforme.
  // void set_time_grid(std::vector<double> time_points);

private:
  rw::models::Gbm model_;
  rw::config::McConfig cfg_;
  bool log_enabled_;
};

} // namespace pricing
} // namespace rw
