#include "rw/calibration/svi.hpp"
#include "rw/pricing/implied_vol.hpp"   // si besoin d’IV pour QC ailleurs
#include <cmath>
#include <algorithm>
#include <limits>

namespace {

inline bool fin(double x){ return std::isfinite(x); }

static void clamp_params(rw::svi::Params& p){
  if (!(fin(p.a)))     p.a = 0.01;
  if (!(fin(p.b)))     p.b = 0.10;
  if (!(fin(p.rho)))   p.rho = 0.0;
  if (!(fin(p.m)))     p.m = 0.0;
  if (!(fin(p.sigma))) p.sigma = 0.20;

  p.b     = std::max(p.b,     1e-8);
  p.sigma = std::max(p.sigma, 1e-6);
  if (p.rho >  0.999) p.rho =  0.999;
  if (p.rho < -0.999) p.rho = -0.999;
}

static double safe_sqrt(double x){
  return std::sqrt(std::max(0.0, x));
}

} // namespace

namespace rw::svi {

double w_svi(double k, const Params& p){
  if (!fin(k)) return std::numeric_limits<double>::quiet_NaN();
  if (!fin(p.a) || !fin(p.b) || !fin(p.rho) || !fin(p.m) || !fin(p.sigma)) return std::numeric_limits<double>::quiet_NaN();
  if (p.b < 0.0 || p.sigma <= 0.0 || std::abs(p.rho) >= 1.0) return std::numeric_limits<double>::quiet_NaN();
  const double z = k - p.m;
  const double s = safe_sqrt(z*z + p.sigma*p.sigma);
  return p.a + p.b * ( p.rho * z + s );
}

svi::FitReport fit_slice(const std::vector<double>& k,
                         const std::vector<double>& w,
                         Params& out,
                         int max_iters,
                         double tol_step)
{
  FitReport rep{};
  rep.n = std::min(k.size(), w.size());
  if (rep.n < 3) { // pas assez de points pour 5 paramètres → fallback “plat”
    out = Params{};
    rep.converged = false;
    rep.rmse_w = std::numeric_limits<double>::quiet_NaN();
    return rep;
  }

  // Init heuristique
  double w_min = std::numeric_limits<double>::infinity();
  double w_max = -std::numeric_limits<double>::infinity();
  double x_med = 0.0;
  {
    std::vector<double> xs = k, ws = w;
    std::sort(xs.begin(), xs.end());
    x_med = xs[xs.size()/2];
    for (std::size_t i=0;i<rep.n;++i){
      if (!fin(w[i]) || !fin(k[i])) continue;
      w_min = std::min(w_min, w[i]);
      w_max = std::max(w_max, w[i]);
    }
    if (!fin(w_min)) w_min = 0.01;
    if (!fin(w_max)) w_max = w_min + 0.05;
  }

  Params p{ std::max(1e-6, w_min*0.8),                      // a
            std::max(1e-4, (w_max - w_min) * 0.5),          // b
            0.0,                                            // rho
            x_med,                                          // m
            0.30 };                                         // sigma

  clamp_params(p);

  auto objective = [&](const Params& pr)->double{
    long double ss = 0.0L;
    for (std::size_t i=0;i<rep.n;++i){
      if (!fin(k[i]) || !fin(w[i])) continue;
      double wi = w_svi(k[i], pr);
      if (!fin(wi)) continue;
      long double r = (long double)wi - (long double)w[i];
      ss += r*r;
    }
    return (double)ss;
  };

  double f_prev = objective(p);
  if (!fin(f_prev)) f_prev = std::numeric_limits<double>::infinity();

  // Levenberg "maison": (J^T J + λ I) Δ = -J^T r, λ adaptatif
  double lambda = 1e-2;
  int it=0; bool ok=false;

  for (; it<max_iters; ++it){
    // Accumulation J^T J et J^T r
    long double A[5][5] = {{0}}; // 5x5
    long double g[5]    = {0,0,0,0,0}; // -J^T r (on remplira avec J^T r, on mettra le signe après)
    std::size_t used = 0;

    for (std::size_t i=0;i<rep.n;++i){
      if (!fin(k[i]) || !fin(w[i])) continue;
      const double z = k[i] - p.m;
      const double s = safe_sqrt(z*z + p.sigma*p.sigma);
      if (s <= 0.0) continue;

      const double wi = p.a + p.b * ( p.rho*z + s );
      if (!fin(wi)) continue;

      const double ri = wi - w[i];

      // dérivées:
      double Ja   = 1.0;
      double Jb   = (p.rho*z + s);
      double Jrho = p.b * z;
      double Jm   = -p.b * ( p.rho + z/s );
      double Jsig = p.b * ( p.sigma / s );

      const double J[5] = {Ja, Jb, Jrho, Jm, Jsig};

      for (int r1=0;r1<5;++r1){
        g[r1] += (long double)J[r1] * (long double)ri; // J^T r
        for (int r2=0;r2<5;++r2){
          A[r1][r2] += (long double)J[r1] * (long double)J[r2]; // J^T J
        }
      }
      ++used;
    }

    if (used < 3) break;

    // Ajoute λ I
    for (int d=0; d<5; ++d) A[d][d] += lambda;

    // Résout (A) Δ = -g
    // Petite Gauss (5x5) — c’est minuscule donc OK
    long double M[5][6]; // [A | -g]
    for (int r1=0;r1<5;++r1){
      for (int r2=0;r2<5;++r2) M[r1][r2] = A[r1][r2];
      M[r1][5] = -g[r1];
    }
    // élimination
    for (int col=0; col<5; ++col){
      // pivot max
      int piv = col;
      long double best = std::fabs(M[piv][col]);
      for (int r=col+1;r<5;++r){ long double v = std::fabs(M[r][col]); if (v>best){best=v; piv=r;} }
      if (best < 1e-18L) { ok=false; goto finish; }
      if (piv!=col) for (int c=col;c<=5;++c) std::swap(M[piv][c], M[col][c]);
      long double pv = M[col][col];
      for (int c=col;c<=5;++c) M[col][c] /= pv;
      for (int r=0;r<5;++r){
        if (r==col) continue;
        long double f = M[r][col];
        for (int c=col;c<=5;++c) M[r][c] -= f * M[col][c];
      }
    }
    long double dlt[5];
    for (int i2=0;i2<5;++i2) dlt[i2] = M[i2][5];

    // étape
    Params cand = p;
    cand.a     += (double)dlt[0];
    cand.b     += (double)dlt[1];
    cand.rho   += (double)dlt[2];
    cand.m     += (double)dlt[3];
    cand.sigma += (double)dlt[4];
    clamp_params(cand);

    const double f_new = objective(cand);
    if (!fin(f_new)) { lambda *= 10.0; continue; }

    const double step_norm = std::sqrt(
      (double)(dlt[0]*dlt[0] + dlt[1]*dlt[1] + dlt[2]*dlt[2] + dlt[3]*dlt[3] + dlt[4]*dlt[4])
    );

    if (f_new < f_prev) {
      p = cand;
      if (step_norm < tol_step || std::fabs(f_prev - f_new) <= 1e-18) { ok=true; f_prev=f_new; break; }
      f_prev = f_new;
      lambda = std::max(1e-9, lambda*0.3);
    } else {
      lambda *= 10.0;
    }
  }

finish:
  out = p;
  rep.converged = ok;
  rep.iters = it+1;

  // RMSE w
  {
    long double ss=0.0L; std::size_t cnt=0;
    for (std::size_t i=0;i<rep.n;++i){
      if (!fin(k[i]) || !fin(w[i])) continue;
      double wi = w_svi(k[i], p);
      if (!fin(wi)) continue;
      long double r = (long double)wi - (long double)w[i];
      ss += r*r; ++cnt;
    }
    rep.rmse_w = (cnt? std::sqrt((double)ss/(double)cnt) : std::numeric_limits<double>::quiet_NaN());
  }

  // RMSE iv (indicatif): on compare sqrt(w/T) vs sqrt(w_svi/T) si l’appelant passe des iv de ref ailleurs.
  rep.rmse_iv = std::numeric_limits<double>::quiet_NaN(); // laissé au CLI si besoin

  return rep;
}

double Surface::iv(double K, double T) const {
  if (!(fin(K) && K>0.0 && fin(T) && T>0.0) || slices.empty()) return std::numeric_limits<double>::quiet_NaN();
  // clamp en T
  if (T <= slices.front().T) {
    const double k = std::log(K / S0);
    const double w = w_svi(k, slices.front().p);
    return fin(w) && w>0.0 ? std::sqrt(w / slices.front().T) : std::numeric_limits<double>::quiet_NaN();
  }
  if (T >= slices.back().T) {
    const double k = std::log(K / S0);
    const double w = w_svi(k, slices.back().p);
    return fin(w) && w>0.0 ? std::sqrt(w / slices.back().T) : std::numeric_limits<double>::quiet_NaN();
  }
  // bracket
  auto it = std::upper_bound(
    slices.begin(), slices.end(), T,
    [](double t, const Slice& s){ return t < s.T; }
  );
  std::size_t j = std::max<std::size_t>(1, it - slices.begin()) - 1;
  const Slice& s1 = slices[j];
  const Slice& s2 = slices[j+1];

  const double k = std::log(K / S0);
  const double w1 = w_svi(k, s1.p);
  const double w2 = w_svi(k, s2.p);
  if (!(fin(w1) && fin(w2))) return std::numeric_limits<double>::quiet_NaN();

  const double a = (T - s1.T) / (s2.T - s1.T);
  const double w = (1.0-a)*w1 + a*w2;   // linéaire en w pour éviter le calendar-arb local
  return w>0.0 ? std::sqrt(w/T) : std::numeric_limits<double>::quiet_NaN();
}

} // namespace rw::svi
