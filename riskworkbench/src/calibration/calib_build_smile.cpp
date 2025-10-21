#include "rw/calibration/calib.hpp"
#include "rw/market/surface.hpp"
#include "rw/smile/smile.hpp"
#include "rw/pricing/implied_vol.hpp"
#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

namespace {
inline bool fin(double x){ return std::isfinite(x); }
}

namespace rw::calib {

rw::smile::SmileSurface build_smile_surface(const rw::market::MarketSurface& mkt,
                                            double tolT)
{
  rw::smile::SmileSurface surf;
  surf.S0 = mkt.S0; surf.r = mkt.r; surf.q = mkt.q;

  // Regroupe par maturité avec tolérance |ΔT| < tolT
  std::map<long long, std::vector<const rw::io::QuoteRow*>> groups;
  const double qTol = std::max(tolT, 1e-12);
  auto keyT = [&](double T)->long long {
    return std::llround(T / qTol);
  };
  for (const auto& r : mkt.rows) {
    if (!(fin(r.K) && fin(r.T)) || r.K<=0.0 || r.T<=0.0) continue;
    groups[keyT(r.T)].push_back(&r);
  }

  for (auto& [Tkq, vec] : groups) {
    if (vec.size() < 2) continue; // pas assez de points
    const double T = Tkq * qTol;

    std::vector<double> xs, ys; xs.reserve(vec.size()); ys.reserve(vec.size());
    for (auto* pr : vec) {
      // IV: prend iv_mid sinon inversion sur mid_price s’il existe
      double iv = std::numeric_limits<double>::quiet_NaN();
      if (fin(pr->iv_mid) && pr->iv_mid>0.0) {
        iv = pr->iv_mid;
      } else if (fin(pr->mid_price)) {
        auto res = rw::bs::implied_vol_cp(mkt.S0, pr->K, mkt.r, mkt.q, pr->T, pr->mid_price, pr->is_call);
        if (res.converged) iv = res.sigma;
      }
      if (!fin(iv) || iv<=0.0) continue;
      xs.push_back(std::log(pr->K / mkt.S0));  // x = ln(K/S0)
      ys.push_back(iv);                        // y = iv
    }

    if (xs.size() < 2) continue;

    // trie par x
    std::vector<std::size_t> idx(xs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto i, auto j){ return xs[i] < xs[j]; });

    // dédoublonne x identiques (moyenne des iv)
    std::vector<double> x2, y2;
    for (std::size_t k=0; k<idx.size(); ) {
      std::size_t j = k+1;
      double sx = xs[idx[k]], sy = ys[idx[k]];
      int cnt=1;
      while (j<idx.size() && std::fabs(xs[idx[j]] - sx) <= 1e-12) {
        sy += ys[idx[j]]; ++cnt; ++j;
      }
      x2.push_back(sx);
      y2.push_back(sy / cnt);
      k = j;
    }
    if (x2.size() < 2) continue;

    rw::smile::Smile1D slice;
    slice.T  = T;
    slice.S0 = mkt.S0;
    slice.x  = std::move(x2);
    slice.y  = std::move(y2);
    slice.build(); // calcule les dérivées PCHIP
    surf.slices.push_back(std::move(slice));
  }

  // ordonne par T
  std::sort(surf.slices.begin(), surf.slices.end(),
            [](const rw::smile::Smile1D& a, const rw::smile::Smile1D& b){ return a.T < b.T; });

  return surf;
}

} // namespace rw::calib
