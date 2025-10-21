#include "rw/market/surface.hpp"
#include "rw/calibration/calib.hpp"
#include "rw/smile/smile.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

static void usage(const char* a){
  std::cerr << "Usage: " << a << " -f market.csv [--check]\n";
}

int main(int argc, char** argv){
  std::string path;
  bool do_check = false;

  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if ((a=="-f"||a=="--file") && i+1<argc) path = argv[++i];
    else if (a=="--check") do_check = true;
    else if (a=="-h"||a=="--help"){ usage(argv[0]); return 0; }
  }
  if (path.empty()) { usage(argv[0]); return 1; }

  rw::market::MarketSurface ms;
  try {
    ms = rw::market::read_market_surface(path);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n"; return 2;
  }

  auto sm = rw::calib::build_smile_surface(ms);

  std::cout << "SmileSurface: slices=" << sm.slices.size()
            << "  S0=" << sm.S0 << " r=" << sm.r << " q=" << sm.q << "\n";

  // Optional self-check: the PCHIP should pass through the nodes
  if (do_check) {
    std::size_t total=0, bad=0;
    for (const auto& sl : sm.slices){
      for (std::size_t i=0;i<sl.x.size();++i){
        const double K = sl.S0 * std::exp(sl.x[i]);
        const double iv_interp = sl.iv_at(K);
        const double iv_node   = sl.y[i];
        if (!(std::isfinite(iv_interp) && std::isfinite(iv_node))) continue;
        ++total;
        if (std::fabs(iv_interp - iv_node) > 1e-8) ++bad;
      }
    }
    std::cout << "Check: pass-through nodes total=" << total
              << " mismatches=" << bad << (bad==0? "  (OK)":"  (FAIL)") << "\n";
  }

  // Sample evaluations: for each slice, print IV at K={0.9S0, S0, 1.1S0}
  for (const auto& sl : sm.slices){
    const double T = sl.T;
    const double K1 = 0.9*sm.S0, K2 = sm.S0, K3 = 1.1*sm.S0;
    const double iv1 = sm.iv_at(T, K1);
    const double iv2 = sm.iv_at(T, K2);
    const double iv3 = sm.iv_at(T, K3);
    std::cout << "T=" << T
              << "  iv(0.9*S0)=" << iv1
              << "  iv(S0)="     << iv2
              << "  iv(1.1*S0)=" << iv3 << "\n";
  }

  return 0;
}
