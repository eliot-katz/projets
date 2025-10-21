#include "rw/calibration/calib.hpp"
#include "rw/market/surface.hpp"
#include "rw/io/market_csv.hpp"
#include "rw/pricing/implied_vol.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace {
inline double df(double r,double T){ return std::exp(-r*T); }
inline double N(double x){ return 0.5 * std::erfc(-x/std::sqrt(2.0)); }
//inline bool std::isfinite(double x){ return std::isstd::isfinite(x); }

static double bs_price(double S, double K, double r, double q, double T, double sigma, bool is_call){
  if (T<=0.0 || sigma<=0.0) {
    const double fwd = S*std::exp((r-q)*T);
    const double intr = is_call ? std::max(0.0, fwd-K) : std::max(0.0, K-fwd);
    return intr * std::exp(-r*T);
  }
  const double v = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T)/v;
  const double d2 = d1 - v;
  const double disc_r = std::exp(-r*T), disc_q = std::exp(-q*T);
  return is_call ? S*disc_q*N(d1) - K*disc_r*N(d2)
                 : K*disc_r*N(-d2) - S*disc_q*N(-d1);
}

static double bs_vega(double S, double K, double r, double q, double T, double sigma){
  if (T<=0.0 || sigma<=0.0) return 0.0;
  const double v = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T)/v;
  const double disc_q = std::exp(-q*T);
  const double nprime = std::exp(-0.5*d1*d1)/std::sqrt(2.0*M_PI);
  return S*disc_q*std::sqrt(T)*nprime;
}

struct Row { double K,T, mid; bool call; };

} // namespace

namespace rw::calib {

CalibReport fit_bs_global_sigma(const rw::market::MarketSurface& mkt,
                                double& sigma_opt,
                                double sigma0,
                                double sigma_min,
                                double sigma_max,
                                int max_iters,
                                double tol_step)
{
  CalibReport rep{};

  // 1) Collecte des quotes "priced" & arbitrage-OK
  std::vector<Row> rows; rows.reserve(mkt.rows.size());
  for (const auto& r : mkt.rows) {
    if (!std::isfinite(r.mid_price) || !std::isfinite(r.K) || !std::isfinite(r.T) || r.K<=0.0 || r.T<=0.0) continue;
    const double lo = r.is_call ? std::max(0.0, mkt.S0*df(mkt.q,r.T) - r.K*df(mkt.r,r.T))
                                : std::max(0.0, r.K*df(mkt.r,r.T) - mkt.S0*df(mkt.q,r.T));
    const double hi = r.is_call ? mkt.S0*df(mkt.q,r.T) : r.K*df(mkt.r,r.T);
    if (r.mid_price < lo-1e-12 || r.mid_price > hi+1e-12) continue;
    rows.push_back({r.K, r.T, r.mid_price, r.is_call});
  }
  rep.n = rows.size();
  if (rows.empty()) { sigma_opt = std::numeric_limits<double>::quiet_NaN(); rep.converged=false; return rep; }

  // 2) Init sigma
  double sigma = sigma0;
  if (!(std::isfinite(sigma) && sigma>0.0)) {
    std::vector<double> ivs;
    for (const auto& r : mkt.rows) if (std::isfinite(r.iv_mid) && r.iv_mid>0.0) ivs.push_back(r.iv_mid);
    if (!ivs.empty()) { std::sort(ivs.begin(), ivs.end()); sigma = ivs[ivs.size()/2]; }
    else sigma = 0.20;
  }
  sigma = std::clamp(sigma, sigma_min, sigma_max);

  auto objective = [&](double s)->double{
    long double ss = 0.0L;
    for (const auto& r: rows) {
      const double p = bs_price(mkt.S0, r.K, mkt.r, mkt.q, r.T, s, r.call);
      const double e = p - r.mid;
      ss += (long double)e*e;
    }
    return (double)ss;
  };

  // 3) Gauss–Newton 1D + backtracking
  int iter=0; bool ok=false; double f_prev = objective(sigma);
  for (; iter<max_iters; ++iter) {
    long double g=0.0L, H=0.0L;
    for (const auto& r : rows) {
      const double p = bs_price(mkt.S0, r.K, mkt.r, mkt.q, r.T, sigma, r.call);
      const double e = p - r.mid;
      const double vega = bs_vega(mkt.S0, r.K, mkt.r, mkt.q, r.T, sigma);
      g += 2.0L * (long double)e * (long double)vega;
      H += 2.0L * (long double)vega * (long double)vega;
    }
    if (H <= 0.0L) break;

    double step = (double)(g / H);
    if (std::fabs(step) < tol_step) { ok=true; break; }

    double sigma_new = std::clamp(sigma - step, sigma_min, sigma_max);
    double f_new = objective(sigma_new);
    int bt=0;
    while (f_new > f_prev && bt<20) {
      step *= 0.5;
      sigma_new = std::clamp(sigma - step, sigma_min, sigma_max);
      f_new = objective(sigma_new);
      ++bt;
    }

    sigma = sigma_new;
    if (std::fabs(step) < tol_step || std::fabs(f_prev - f_new) <= 1e-18) { ok=true; f_prev=f_new; break; }
    f_prev = f_new;
  }

  sigma_opt = sigma;
  rep.converged = ok;
  rep.iters = iter+1;

  // 4) RMSE prix
  {
    long double ss=0.0L;
    for (const auto& r: rows) {
      const double p = bs_price(mkt.S0, r.K, mkt.r, mkt.q, r.T, sigma_opt, r.call);
      const double e = p - r.mid;
      ss += (long double)e*e;
    }
    rep.rmse_price = std::sqrt((double)ss / (double)rows.size());
  }

  // 5) RMSE IV (utilise iv_mid si dispo, sinon inversion BS)
  {
    std::size_t cnt = 0; long double ss=0.0L;
    for (const auto& r : mkt.rows) {
      if (!(std::isfinite(r.K) && std::isfinite(r.T) && r.K>0.0 && r.T>0.0)) continue;
      double iv_i = std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(r.iv_mid) && r.iv_mid>0.0) iv_i = r.iv_mid;
      else if (std::isfinite(r.mid_price)) {
        auto res = rw::bs::implied_vol_cp(mkt.S0, r.K, mkt.r, mkt.q, r.T, r.mid_price, r.is_call);
        if (res.converged) iv_i = res.sigma;
      }
      if (!std::isfinite(iv_i)) continue;
      const double e = iv_i - sigma_opt;
      ss += (long double)e*e;
      ++cnt;
    }
    rep.rmse_iv = (cnt>0) ? std::sqrt((double)ss / (double)cnt)
                          : std::numeric_limits<double>::quiet_NaN();
  }

  return rep;
}

} // namespace rw::calib
