#include "rw/market/surface.hpp"
#include "rw/calibration/calib.hpp"
#include <iostream>
#include <string>
#include <cstdlib>

static void usage(const char* a){
  std::cerr <<
  "Usage: " << a << " -f market.csv [--sigma0 0.2] [--min 1e-6] [--max 5] [--maxit 100]\n";
}

int main(int argc, char** argv){
  std::string path;
  double sigma0 = 0.0;  // 0 => auto init
  double smin = 1e-6, smax=5.0;
  int maxit = 100;

  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if ((a=="-f"||a=="--file") && i+1<argc) path = argv[++i];
    else if (a=="--sigma0" && i+1<argc) sigma0 = std::atof(argv[++i]);
    else if (a=="--min" && i+1<argc) smin = std::atof(argv[++i]);
    else if (a=="--max" && i+1<argc) smax = std::atof(argv[++i]);
    else if (a=="--maxit" && i+1<argc) maxit = std::atoi(argv[++i]);
    else if (a=="-h"||a=="--help") { usage(argv[0]); return 0; }
  }
  if (path.empty()) { usage(argv[0]); return 1; }

  rw::market::MarketSurface ms = rw::market::read_market_surface(path);

  double sigma_opt = 0.0;
  const auto rep = rw::calib::fit_bs_global_sigma(ms, sigma_opt, sigma0, smin, smax, maxit);

  std::cout.setf(std::ios::fixed); std::cout.precision(8);
  std::cout << "Global sigma: " << sigma_opt << (rep.converged? " (converged)":" (not converged)") << "\n"
            << "RMSE_price: " << rep.rmse_price << "\n"
            << "RMSE_iv: " << rep.rmse_iv << "\n"
            << "points_used: " << rep.n << "  iters: " << rep.iters << "\n";

  return rep.converged ? 0 : 2;
}
