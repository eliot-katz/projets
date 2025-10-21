#include "rw/market/iv_surface.hpp"
#include "rw/io/market_csv.hpp"
#include "rw/pricing/implied_vol.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

using std::size_t;

namespace {
inline bool is_finite(double x){ return std::isfinite(x); }
inline double df(double r, double T){ return std::exp(-r*T); }

// (Optionnel, utile si tu veux QC plus tard)
static double bs_call(double S, double K, double r, double q, double T, double sigma){
  if (T<=0.0 || sigma<=0.0) return std::max(0.0, S*df(q,T) - K*df(r,T));
  double v = sigma * std::sqrt(T);
  double m = std::log(S/K) + (r-q)*T;
  double d1 = (m + 0.5*sigma*sigma*T)/v;
  double d2 = d1 - v;
  auto N = [](double x){ return 0.5*std::erfc(-x/std::sqrt(2.0)); };
  return S*df(q,T)*N(d1) - K*df(r,T)*N(d2);
}
} // namespace

namespace rw::market {

// ===== IvSurface =====

IvSurface::IvSurface(double S0, double r, double q)
  : S0_(S0), r_(r), q_(q) {}

void IvSurface::add_quote(double K, double T, bool is_call,
                          std::optional<double> mid_price,
                          std::optional<double> iv_mid)
{
  Pending p{};
  p.K = K; p.T = T; p.is_call = is_call;
  p.has_price = mid_price.has_value() && is_finite(*mid_price);
  p.price     = p.has_price ? *mid_price : std::numeric_limits<double>::quiet_NaN();
  p.has_iv    = iv_mid.has_value() && is_finite(*iv_mid);
  p.iv        = p.has_iv ? *iv_mid : std::numeric_limits<double>::quiet_NaN();
  pending_.push_back(p);
}

void IvSurface::add_quotes_from_rows(const std::vector<rw::io::QuoteRow>& rows) {
  pending_.reserve(pending_.size() + rows.size());
  for (const auto& r : rows) {
    std::optional<double> mp = std::isfinite(r.mid_price) ? std::optional<double>(r.mid_price) : std::nullopt;
    std::optional<double> iv = std::isfinite(r.iv_mid)    ? std::optional<double>(r.iv_mid)    : std::nullopt;
    add_quote(r.K, r.T, r.is_call, mp, iv);
  }
}

// ===== PCHIP (Fritsch–Carlson) =====
void IvSurface::Slice::build() {
  d.clear();
  const size_t n = x.size();
  if (n == 0) { built = true; return; }
  if (n == 1) { d.assign(1, 0.0); built = true; return; }

  std::vector<double> h(n-1), delta(n-1);
  for (size_t i=0;i+1<n;++i) {
    h[i] = x[i+1]-x[i];
    delta[i] = (w[i+1]-w[i]) / h[i];
  }

  d.assign(n, 0.0);

  // endpoints
  if (n >= 2) {
    if (n==2) {
      d[0] = d[1] = delta[0];
    } else {
      double d0 = ((2*h[0]+h[1])*delta[0] - h[0]*delta[1]) / (h[0]+h[1]);
      if (d0*delta[0] <= 0) d0 = 0.0;
      else if (std::fabs(d0) > 2*std::fabs(delta[0])) d0 = 2*delta[0];
      d[0] = d0;

      double dn = ((2*h[n-2]+h[n-3])*delta[n-2] - h[n-2]*delta[n-3]) / (h[n-3]+h[n-2]);
      if (dn*delta[n-2] <= 0) dn = 0.0;
      else if (std::fabs(dn) > 2*std::fabs(delta[n-2])) dn = 2*delta[n-2];
      d[n-1] = dn;

      for (size_t i=1;i+1<n;++i) {
        if (delta[i-1]*delta[i] <= 0.0) { d[i] = 0.0; }
        else {
          double w1 = 2*h[i] + h[i-1];
          double w2 = h[i] + 2*h[i-1];
          d[i] = (w1+w2) / (w1/delta[i-1] + w2/delta[i]);
        }
      }
    }
  }
  built = true;
}

double IvSurface::Slice::eval(double xq) const {
  const size_t n = x.size();
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  if (n == 1) return w[0];

  // clamp extrapolation (constant aux bords)
  if (xq <= x.front()) return w.front();
  if (xq >= x.back())  return w.back();

  auto it = std::upper_bound(x.begin(), x.end(), xq);
  size_t i = std::max<size_t>(1, it - x.begin()) - 1;

  double h = x[i+1]-x[i];
  double t = (xq - x[i]) / h;
  double h00 = (2*t*t*t - 3*t*t + 1);
  double h10 = (t*t*t - 2*t*t + t);
  double h01 = (-2*t*t*t + 3*t*t);
  double h11 = (t*t*t - t*t);

  return h00*w[i] + h10*h*d[i] + h01*w[i+1] + h11*h*d[i+1];
}

double IvSurface::clamp(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

void IvSurface::build() {
  kept_ = skipped_ = 0;

  // 1) Convertit toutes les quotes en (T, x, w), en tentant price->IV si besoin
  struct Node { double T; double x; double w; };
  std::vector<Node> nodes; nodes.reserve(pending_.size());

  for (const auto& p : pending_) {
    if (!(is_finite(p.K) && p.K > 0.0 && is_finite(p.T) && p.T > 0.0)) { ++skipped_; continue; }

    double iv = std::numeric_limits<double>::quiet_NaN();
    if (p.has_iv) {
      iv = p.iv;
    } else if (p.has_price) {
      auto res = rw::bs::implied_vol_cp(S0_, p.K, r_, q_, p.T, p.price, p.is_call);
      if (res.converged) iv = res.sigma;
    }

    if (!is_finite(iv) || iv <= 0.0) { ++skipped_; continue; }
    double w = iv*iv * p.T;
    if (!is_finite(w) || w <= 0.0) { ++skipped_; continue; }

    nodes.push_back(Node{p.T, std::log(p.K / S0_), w});
    ++kept_;
  }

  pending_.clear();

  if (nodes.empty()) { slices_.clear(); return; }

  // 2) Regroupe par maturité (clé T)
  std::map<double, std::vector<Node>> byT;
  for (auto& n : nodes) byT[n.T].push_back(n);

  // 3) Crée les slices triées par x
  slices_.clear();
  slices_.reserve(byT.size());
  for (auto& kv : byT) {
    Slice s; s.T = kv.first;
    auto& vec = kv.second;
    std::sort(vec.begin(), vec.end(), [](const Node& a, const Node& b){ return a.x < b.x; });

    // dédoublonne x identiques (moyenne de w)
    std::vector<double> xs; xs.reserve(vec.size());
    std::vector<double> ws; ws.reserve(vec.size());
    for (size_t i=0;i<vec.size();){
      double x = vec[i].x;
      double w = vec[i].w;
      size_t j = i+1;
      while (j < vec.size() && std::fabs(vec[j].x - x) < 1e-12) { w = 0.5*(w + vec[j].w); ++j; }
      xs.push_back(x); ws.push_back(w);
      i = j;
    }
    s.x = std::move(xs);
    s.w = std::move(ws);
    s.build();
    slices_.push_back(std::move(s));
  }

  // 4) Trie les slices par T
  std::sort(slices_.begin(), slices_.end(), [](const Slice& a, const Slice& b){ return a.T < b.T; });
}

std::size_t IvSurface::maturities() const { return slices_.size(); }

// Interpolation en T : linéaire sur w entre T-brackets (clamp si hors domaine)
double IvSurface::iv(double K, double T) const {
  if (slices_.empty() || !(is_finite(K) && K>0.0) || !(is_finite(T) && T>0.0))
    return std::numeric_limits<double>::quiet_NaN();

  const double xq = std::log(K / S0_);

  // T hors domaine : clamp aux slices extrêmes
  if (T <= slices_.front().T) {
    double w = slices_.front().eval(xq);
    return std::sqrt(std::max(0.0, w / std::max(T, 1e-16)));
  }
  if (T >= slices_.back().T) {
    double w = slices_.back().eval(xq);
    return std::sqrt(std::max(0.0, w / T));
  }

  // cherche bracket en T
  auto it = std::upper_bound(slices_.begin(), slices_.end(), T,
    [](double t, const Slice& s){ return t < s.T; });
  size_t j = std::max<size_t>(1, it - slices_.begin()) - 1;
  const Slice& s1 = slices_[j];
  const Slice& s2 = slices_[j+1];

  const double w1 = s1.eval(xq);
  const double w2 = s2.eval(xq);
  const double alpha = (T - s1.T) / (s2.T - s1.T);
  const double w = (1.0 - alpha) * w1 + alpha * w2;

  return std::sqrt(std::max(0.0, w / T));
}

} // namespace rw::market
