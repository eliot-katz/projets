#include "rw/smile/smile.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

inline bool fin(double x){ return std::isfinite(x); }

// PCHIP: calcule les pentes d[i] (Fritsch–Carlson) à partir de x[], y[] strictement croissants en x.
static void pchip_slopes(const std::vector<double>& x,
                         const std::vector<double>& y,
                         std::vector<double>& d)
{
  const std::size_t n = x.size();
  d.assign(n, 0.0);
  if (n < 2) return;

  // h_i, delta_i
  std::vector<double> h(n-1), delta(n-1);
  for (std::size_t i=0;i+1<n;++i){
    h[i] = x[i+1] - x[i];
    delta[i] = (y[i+1] - y[i]) / h[i];
  }

  // Extrémités: dérivées une–sided “shape-preserving”
  d[0] = delta[0];
  d[n-1] = delta[n-2];

  // Intérieur
  for (std::size_t i=1;i+1<n;++i){
    if (delta[i-1] * delta[i] <= 0.0) {
      d[i] = 0.0;
    } else {
      const double w1 = 2.0*h[i] + h[i-1];
      const double w2 = h[i] + 2.0*h[i-1];
      d[i] = (w1 + w2 > 0.0) ? (w1 + w2) / (w1/delta[i-1] + w2/delta[i]) : 0.0;
    }
  }

  // Safeguard (Fritsch–Butland): clamp quand pente trop “raide” par rapport à delta
  for (std::size_t i=0;i+1<n;++i){
    if (delta[i] == 0.0) { d[i] = 0.0; d[i+1] = 0.0; }
    else {
      const double a = d[i]   / delta[i];
      const double b = d[i+1] / delta[i];
      const double s = a*a + b*b;
      if (s > 9.0) {
        const double t = 3.0 / std::sqrt(s);
        d[i]   = t * a * delta[i];
        d[i+1] = t * b * delta[i];
      }
    }
  }
}

// Évalue l’Hermite cubique sur [x0,x1] à t∈[0,1] (y0,y1 pentes m0,m1).
inline double hermite(double t, double y0, double y1, double m0, double m1, double h){
  const double t2 = t*t, t3 = t2*t;
  const double h00 =  2*t3 - 3*t2 + 1;
  const double h10 =    t3 - 2*t2 + t;
  const double h01 = -2*t3 + 3*t2;
  const double h11 =    t3 -   t2;
  return h00*y0 + h10*h*m0 + h01*y1 + h11*h*m1;
}

} // namespace

namespace rw::smile {

void Smile1D::build(){
  d.clear();
  if (x.size() < 2 || x.size()!=y.size()) return;
  // suppose x strictement croissante (déjà nettoyé en amont)
  pchip_slopes(x, y, d);
}

double Smile1D::iv_at_x(double xq) const {
  const std::size_t n = x.size();
  if (n < 2 || n != y.size() || d.size()!=n || !fin(xq))
    return std::numeric_limits<double>::quiet_NaN();

  // bornes raisonnables pour l'IV (éviter négatifs/absurdités en extrapolation)
  constexpr double IV_MIN = 1e-8;
  constexpr double IV_MAX = 5.0;

  auto clamp_iv = [&](double iv){
    if (!fin(iv)) return std::numeric_limits<double>::quiet_NaN();
    if (iv < IV_MIN) iv = IV_MIN;
    if (iv > IV_MAX) iv = IV_MAX;
    return iv;
  };

  // En-dehors: extrapolation PLATE (valeur du bord), beaucoup plus stable
  if (xq <= x.front()) {
    return clamp_iv(y.front());
  }
  if (xq >= x.back()) {
    return clamp_iv(y.back());
  }

  // Cherche l’intervalle [i, i+1]
  auto it = std::upper_bound(x.begin(), x.end(), xq);
  std::size_t i = std::max<std::size_t>(1, it - x.begin()) - 1;

  const double h = x[i+1] - x[i];
  const double t = (xq - x[i]) / h;
  const double iv = hermite(t, y[i], y[i+1], d[i], d[i+1], h);

  return clamp_iv(iv);
}


double Smile1D::iv_at(double K) const {
  if (!(fin(K) && K>0.0 && fin(S0) && S0>0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double xq = std::log(K / S0);
  return iv_at_x(xq);
}

double SmileSurface::iv_at(double T, double K) const {
  if (!(fin(T) && T>0.0 && fin(K) && K>0.0) || slices.empty()) return std::numeric_limits<double>::quiet_NaN();

  // cas bord: clamp en T
  if (T <= slices.front().T) {
    const double iv = slices.front().iv_at(K);
    return fin(iv) ? iv : std::numeric_limits<double>::quiet_NaN();
  }
  if (T >= slices.back().T) {
    const double iv = slices.back().iv_at(K);
    return fin(iv) ? iv : std::numeric_limits<double>::quiet_NaN();
  }

  // bracket [j, j+1] tel que T∈(T_j, T_{j+1}]
  auto it = std::upper_bound(
    slices.begin(), slices.end(), T,
    [](double t, const Smile1D& s){ return t < s.T; }
  );
  std::size_t j = std::max<std::size_t>(1, it - slices.begin()) - 1;
  const Smile1D& s1 = slices[j];
  const Smile1D& s2 = slices[j+1];

  const double iv1 = s1.iv_at(K);
  const double iv2 = s2.iv_at(K);
  if (!(fin(iv1) && fin(iv2))) return std::numeric_limits<double>::quiet_NaN();

  // interp linéaire en w = σ² T pour limiter le calendar arb local
  const double w1 = iv1*iv1 * s1.T;
  const double w2 = iv2*iv2 * s2.T;
  const double a  = (T - s1.T) / (s2.T - s1.T);
  const double w  = (1.0 - a)*w1 + a*w2;
  return (w>0.0) ? std::sqrt(w / T) : std::numeric_limits<double>::quiet_NaN();
}

} // namespace rw::smile
