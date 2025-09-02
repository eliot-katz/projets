#pragma once
/**
 * @file analytic_bs.hpp
 * @brief Formules fermées Black–Scholes (référence & tests).
 *
 * # Modèle (mesure Q)
 * dS_t / S_t = (r - q) dt + sigma dW_t,  avec r,q en décimal (peuvent être négatifs).
 *
 * # Notations
 * d1 = [ ln(S0/K) + (r - q + 0.5*sigma^2) T ] / (sigma * sqrt(T))
 * d2 = d1 - sigma * sqrt(T)
 *
 * Prix (valeurs au temps 0) :
 *   Call = S0 * e^{-qT} * N(d1) - K * e^{-rT} * N(d2)
 *   Put  = K * e^{-rT} * N(-d2) - S0 * e^{-qT} * N(-d1)
 *
 * # Stabilité numérique
 * - Utiliser des branches pour T=0 ou sigma=0 afin d’éviter 0/0 et inf/NaN.
 *   * Si T == 0 : prix = payoff (pas d’actualisation).
 *   * Si sigma == 0 : S_T = S0 * exp((r-q)T) déterministe ⇒ prix = e^{-rT} * payoff(S_T).
 * - Calculer d1/d2 de façon stable (regrouper ln(S0/K) et drift * T).
 *
 * # Unités
 * - T en années fractionnelles (ex : 0.5 = 6 mois).
 * - Taux continus annualisés en décimal.
 *
 * # Tests
 * - Cas “or” (exemples de référence) et tolérances numériques.
 * - Vérifier la parité put–call : C - P == S0 e^{-qT} - K e^{-rT}.
 */

#include <utility> // autorisé (utile si besoin de couples plus tard)

namespace rw {
namespace pricing {

/**
 * @brief Prix Black–Scholes d’un call européen.
 * @param S0    Spot initial (>0)
 * @param K     Strike (>0)
 * @param r     Taux sans risque (décimal, peut être < 0)
 * @param q     Taux de dividende (décimal, peut être < 0)
 * @param sigma Volatilité (>= 0)
 * @param T     Maturité en années (>= 0)
 * @return Prix au temps 0 (double)
 *
 * Cas limites gérés (recommandé dans l’implémentation .cpp) :
 * - T == 0 : max(S0 - K, 0)
 * - sigma == 0 : e^{-rT} * max(S0 * e^{(r-q)T} - K, 0)
 */
double price_call_bs(double S0, double K, double r, double q, double sigma, double T) noexcept;

/**
 * @brief Prix Black–Scholes d’un put européen.
 * @param S0    Spot initial (>0)
 * @param K     Strike (>0)
 * @param r     Taux sans risque (décimal, peut être < 0)
 * @param q     Taux de dividende (décimal, peut être < 0)
 * @param sigma Volatilité (>= 0)
 * @param T     Maturité en années (>= 0)
 * @return Prix au temps 0 (double)
 *
 * Cas limites gérés (recommandé dans l’implémentation .cpp) :
 * - T == 0 : max(K - S0, 0)
 * - sigma == 0 : e^{-rT} * max(K - S0 * e^{(r-q)T}, 0)
 */
double price_put_bs(double S0, double K, double r, double q, double sigma, double T) noexcept;

/**
 * @brief Écart de parité put–call.
 * @param call Prix du call
 * @param put  Prix du put
 * @param S0   Spot initial
 * @param K    Strike
 * @param r    Taux sans risque (décimal)
 * @param q    Taux de dividende (décimal)
 * @param T    Maturité (années)
 * @return  gap = call - put - ( S0 * e^{-qT} - K * e^{-rT} )
 *         (doit être ≈ 0 en BS sans frictions)
 */
double put_call_parity_gap(double call, double put,
                           double S0, double K, double r, double q, double T) noexcept;

} // namespace pricing
} // namespace rw
