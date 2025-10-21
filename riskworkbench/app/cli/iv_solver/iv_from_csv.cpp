// app/cli/iv_solver/iv_from_csv.cpp
#include "rw/io/market_csv.hpp"
#include "rw/pricing/implied_vol.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

static inline double df(double r, double T){ return std::exp(-r*T); }

static bool call_in_bounds(double S0,double K,double r,double q,double T,double C){
  const double lo = std::max(0.0, S0*df(q,T) - K*df(r,T));
  const double hi = S0*df(q,T);
  return C >= lo && C <= hi;
}
static bool put_in_bounds(double S0,double K,double r,double q,double T,double P){
  const double lo = std::max(0.0, K*df(r,T) - S0*df(q,T));
  const double hi = K*df(r,T);
  return P >= lo && P <= hi;
}

int main(int argc, char** argv) {
  std::string path = (argc>1 ? argv[1] : "data/market_samples/sample_market.csv");

  std::size_t ignored=0; std::vector<std::string> warns;
  auto rows = rw::io::read_market_csv(path, &ignored, &warns);

  double S0=NAN,r=NAN,q=NAN;
  if (auto m = rw::io::read_market_meta_csv(path)) { S0=m->S0; r=m->r; q=m->q; }
  if (!std::isfinite(S0)) S0=5000.0;
  if (!std::isfinite(r))  r =0.02;
  if (!std::isfinite(q))  q =0.01;

  int total=0, priced=0, priced_ok=0, conv=0, calls=0, puts=0, priced_bad=0;

  for (const auto& row : rows) {
    ++total;
    if (!std::isfinite(row.mid_price)) continue; // ici on teste seulement les lignes prix
    ++priced;
    row.is_call ? ++calls : ++puts;

    const bool ok = row.is_call
      ? call_in_bounds(S0,row.K,r,q,row.T,row.mid_price)
      : put_in_bounds (S0,row.K,r,q,row.T,row.mid_price);

    if (!ok) { ++priced_bad; continue; } // on ne tente pas l'inversion
    ++priced_ok;

    auto res = rw::bs::implied_vol_cp(S0,row.K,r,q,row.T,row.mid_price,row.is_call);
    if (res.converged) ++conv;
  }

  std::cout << "Rows=" << total
            << " priced=" << priced
            << " priced_ok=" << priced_ok
            << " priced_bad=" << priced_bad
            << " converged=" << conv
            << " (" << (priced_ok? (100.0*conv/priced_ok):0.0) << "% of arbitrage-OK)"
            << "  calls=" << calls << " puts=" << puts
            << "\n";

  for (auto& w : warns) std::cerr << "[warn] " << w << "\n";
  return 0;
}
