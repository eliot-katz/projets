#include "rw/pricing/implied_vol.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace {

// N(x) et phi(x)
inline double norm_cdf(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}
inline double norm_pdf(double x) {
  static constexpr double INV_SQRT_2PI = 0.39894228040143267793994605993438; // 1/sqrt(2π)
  return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

// Prix BS call/put + vega(call)
inline double bs_call(double S, double K, double r, double q, double T, double sigma) {
  if (T <= 0.0 || sigma <= 0.0) {
    const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
    return std::max(0.0, S * df_q - K * df_r);
  }
  const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
  const double volSqrtT = sigma * std::sqrt(T);
  const double m = std::log(S / K) + (r - q) * T;
  const double d1 = (m + 0.5 * sigma * sigma * T) / volSqrtT;
  const double d2 = d1 - volSqrtT;
  return S * df_q * norm_cdf(d1) - K * df_r * norm_cdf(d2);
}
inline double bs_put(double S, double K, double r, double q, double T, double sigma) {
  if (T <= 0.0 || sigma <= 0.0) {
    const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
    return std::max(0.0, K * df_r - S * df_q);
  }
  const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
  const double volSqrtT = sigma * std::sqrt(T);
  const double m = std::log(S / K) + (r - q) * T;
  const double d1 = (m + 0.5 * sigma * sigma * T) / volSqrtT;
  const double d2 = d1 - volSqrtT;
  return K * df_r * norm_cdf(-d2) - S * df_q * norm_cdf(-d1);
}
inline double bs_vega_call(double S, double K, double r, double q, double T, double sigma) {
  if (T <= 0.0 || sigma <= 0.0) return 0.0;
  const double df_q = std::exp(-q * T);
  const double volSqrtT = sigma * std::sqrt(T);
  const double m = std::log(S / K) + (r - q) * T;
  const double d1 = (m + 0.5 * sigma * sigma * T) / volSqrtT;
  return S * df_q * std::sqrt(T) * norm_pdf(d1); // dC/dσ
}

// Parité put-call : C = P + S*e^{-qT} - K*e^{-rT}
inline double put_to_call(double put, double S, double K, double r, double q, double T) {
  const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
  return put + S * df_q - K * df_r;
}

// Bornes de non-arbitrage pour CALL
inline double call_lower_bound(double S, double K, double r, double q, double T) {
  const double df_r = std::exp(-r * T), df_q = std::exp(-q * T);
  return std::max(0.0, S * df_q - K * df_r);
}
inline double call_upper_bound(double S, double /*K*/, double /*r*/, double q, double T) {
  const double df_q = std::exp(-q * T);
  return S * df_q; // σ→∞
}

struct SolveCfg {
  double sigma_lo = 1e-6;
  double sigma_hi = 5.0;
  int    newton_max = 50;
  int    bisect_max = 100;
  double tol_price  = 1e-12;  // tol. absolue sur f(σ) en unités de prix
  double tol_sigma  = 1e-10;  // tol. absolue sur σ
};

} // namespace

namespace rw::bs {

IvResult implied_vol_cp(double S0, double K, double r, double q, double T, double price, bool is_call) {
  // Transforme en prix CALL si besoin
  const double call_price = is_call ? price : put_to_call(price, S0, K, r, q, T);

  SolveCfg cfg{};
  // Vérifs de bornes (non-arbitrage)
  const double c_lo = call_lower_bound(S0, K, r, q, T);
  const double c_hi = call_upper_bound(S0, K, r, q, T);

  // Si le prix est hors bornes, on renvoie la borne correspondante (non convergé)
  if (!(call_price >= c_lo && call_price <= c_hi)) {
    IvResult out{};
    out.sigma     = (call_price < c_lo ? cfg.sigma_lo : cfg.sigma_hi);
    out.iters     = 0;
    out.converged = false;
    return out;
  }

  // Bracketing initial
  double a = cfg.sigma_lo, b = cfg.sigma_hi;
  auto f = [&](double s) { return bs_call(S0, K, r, q, T, s) - call_price; };
  double fa = f(a), fb = f(b);

  // En théorie monotone croissant → fa <= 0 <= fb
  if (fa > 0.0) { fa = f(a = std::max(1e-12, 0.5 * a)); } // tente plus petit
  if (fb < 0.0) { fb = f(b = b * 2.0); }                  // tente plus grand (rarement utile)
  // Si toujours pas un vrai bracket, on continue quand même avec safeguarded Newton (f(a)*f(b) devrait être <=0)

  // Guess initial type Brenner-Subrahmanyam (approx ATM)
  const double df_q = std::exp(-q * T);
  double sigma = std::clamp(std::sqrt(2.0 * M_PI / std::max(T, 1e-16)) * (call_price / std::max(S0 * df_q, 1e-16)),
                            cfg.sigma_lo, cfg.sigma_hi);

  // Safeguarded Newton: si step hors [a,b] ou vega trop petit → bisection
  int it = 0;
  bool conv = false;
  for (; it < cfg.newton_max; ++it) {
    const double price_sigma = bs_call(S0, K, r, q, T, sigma);
    const double diff = price_sigma - call_price;
    if (std::fabs(diff) <= cfg.tol_price * std::max(1.0, call_price)) { conv = true; break; }

    const double vega = bs_vega_call(S0, K, r, q, T, sigma);
    // Met à jour le bracket
    if (diff < 0.0) { a = sigma; fa = diff; } else { b = sigma; fb = diff; }

    // Essaye un pas Newton si vega raisonnable
    double sigma_newton = sigma;
    if (vega > 1e-14) {
      sigma_newton = sigma - diff / vega;
    } else {
      sigma_newton = std::numeric_limits<double>::quiet_NaN(); // forcera bisection
    }

    // Choix du pas
    if (std::isfinite(sigma_newton) && sigma_newton > a && sigma_newton < b) {
      // Safeguard : si amélioration insuffisante, fallback bisection
      const double mid = 0.5 * (a + b);
      if (std::fabs(sigma_newton - sigma) < 0.1 * std::fabs(mid - sigma)) {
        sigma = sigma_newton;
      } else {
        sigma = mid;
      }
    } else {
      // Bisection
      sigma = 0.5 * (a + b);
    }

    if (std::fabs(b - a) < cfg.tol_sigma) { conv = true; break; }
  }

  // Si non convergé en Newton → bisection pure
  if (!conv) {
    for (int k = 0; k < cfg.bisect_max; ++k, ++it) {
      sigma = 0.5 * (a + b);
      const double val = f(sigma);
      if (std::fabs(val) <= cfg.tol_price * std::max(1.0, call_price)) { conv = true; break; }
      if (val < 0.0) { a = sigma; } else { b = sigma; }
      if (std::fabs(b - a) < cfg.tol_sigma) { conv = true; break; }
    }
  }

  IvResult out{};
  out.sigma     = std::clamp(sigma, cfg.sigma_lo, cfg.sigma_hi);
  out.iters     = it;
  out.converged = conv;
  return out;
}

IvResult implied_vol(double S0, double K, double r, double q, double T, double price) {
  // Par défaut, interprète "price" comme un **CALL**
  return implied_vol_cp(S0, K, r, q, T, price, /*is_call=*/true);
}

} // namespace rw::bs
