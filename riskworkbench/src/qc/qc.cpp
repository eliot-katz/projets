#include "rw/qc/qc.hpp"
#include "rw/market/surface.hpp"
#include "rw/smile/smile.hpp"
#include "rw/pricing/implied_vol.hpp"
#include <cmath>
#include <limits>

namespace {
inline bool fin(double x){ return std::isfinite(x); }
inline double N(double x){ return 0.5 * std::erfc(-x * M_SQRT1_2); }

inline double bs_price_cp(double S0, double K, double r, double q, double T,
                          double sigma, bool is_call)
{
  if (!(fin(S0)&&fin(K)&&fin(r)&&fin(q)&&fin(T)&&fin(sigma)) || S0<=0.0 || K<=0.0 || T<0.0)
    return std::numeric_limits<double>::quiet_NaN();

  if (T==0.0 || sigma<=0.0) {
    const double df_r = std::exp(-r*T);
    const double df_q = std::exp(-q*T);
    const double fwd  = S0*df_q - K*df_r;
    return is_call ? std::max(0.0, fwd) : std::max(0.0, -fwd);
  }

  const double df_r = std::exp(-r*T);
  const double df_q = std::exp(-q*T);
  const double vs   = sigma * std::sqrt(T);
  const double m    = std::log(S0/K) + (r - q + 0.5*sigma*sigma) * T;
  const double d1   = m / vs;
  const double d2   = d1 - vs;

  return is_call ? (S0*df_q*N(d1) - K*df_r*N(d2))
                 : (K*df_r*N(-d2) - S0*df_q*N(-d1));
}
} // namespace

namespace rw::qc {

std::vector<ErrorRow>
compare_fit(const rw::market::MarketSurface& mkt,
            const rw::smile::SmileSurface& smile)
{
  std::vector<ErrorRow> out;
  out.reserve(mkt.rows.size());

  for (const auto& r : mkt.rows) {
    if (!(fin(r.K) && fin(r.T)) || r.K<=0.0 || r.T<=0.0) continue;

    // 1) IV marché
    double iv_mkt = std::numeric_limits<double>::quiet_NaN();
    if (fin(r.iv_mid) && r.iv_mid>0.0) {
      iv_mkt = r.iv_mid;
    } else if (fin(r.mid_price)) {
      auto ivres = rw::bs::implied_vol_cp(mkt.S0, r.K, mkt.r, mkt.q, r.T, r.mid_price, r.is_call);
      if (ivres.converged && ivres.sigma>0.0) iv_mkt = ivres.sigma;
    }
    if (!(fin(iv_mkt) && iv_mkt>0.0)) continue;

    // 2) Prix marché
    double price_mkt = std::numeric_limits<double>::quiet_NaN();
    if (fin(r.mid_price)) {
      price_mkt = r.mid_price;
    } else {
      price_mkt = bs_price_cp(mkt.S0, r.K, mkt.r, mkt.q, r.T, iv_mkt, r.is_call);
    }
    if (!fin(price_mkt)) continue;

    // 3) IV fit via smile(T,K) + 4) Prix fit
    const double iv_fit = smile.iv_at(r.T, r.K);
    if (!(fin(iv_fit) && iv_fit>0.0)) continue;

    const double price_fit = bs_price_cp(mkt.S0, r.K, mkt.r, mkt.q, r.T, iv_fit, r.is_call);
    if (!fin(price_fit)) continue;

    ErrorRow e;
    e.K = r.K; e.T = r.T; e.is_call = r.is_call;
    e.price_mkt = price_mkt;
    e.price_fit = price_fit;
    e.err_price = price_fit - price_mkt;
    e.iv_mkt    = iv_mkt;
    e.iv_fit    = iv_fit;
    e.err_iv    = iv_fit - iv_mkt;
    out.push_back(e);
  }

  return out;
}

} // namespace rw::qc
