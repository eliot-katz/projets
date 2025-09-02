#include <rw/models/gbm.hpp>
#include <rw/core/normal_rng.hpp>
#include <rw/core/stats.hpp>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

int main(int argc, char** argv) {
  if (argc != 8) {
    std::cerr << "Usage: " << argv[0] << " S0 r q sigma T seed N\n";
    return 1;
  }
  double S0   = std::stod(argv[1]);
  double r    = std::stod(argv[2]);
  double q    = std::stod(argv[3]);
  double sig  = std::stod(argv[4]);
  double T    = std::stod(argv[5]);
  unsigned long long seed = std::stoull(argv[6]);
  std::size_t N = static_cast<std::size_t>(std::stoull(argv[7]));

  rw::models::Gbm model({r,q,sig});
  rw::core::NormalRng rng(seed);

  // Accumulateurs pour S_T (pas d'actualisation)
  rw::core::RunningStats acc;
  for (std::size_t i=0; i<N; ++i) {
    double ST = model.sample_ST(S0, T, rng);
    acc.add(ST);
  }

  double mean_emp = acc.mean();
  // variance d'échantillon acc.variance(); on préfère une variance "pop" ici
  // mais sur N grand ça change peu :
  double var_emp = acc.variance();

  double mean_th = S0 * std::exp((r - q) * T);
  double var_th  = mean_th * mean_th * (std::exp(sig*sig*T) - 1.0);

  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(6);
  std::cout << "mean_emp=" << mean_emp << " mean_th=" << mean_th
            << " rel_err=" << std::abs(mean_emp-mean_th)/mean_th << "\n";
  std::cout << "var_emp="  << var_emp  << " var_th="  << var_th
            << " rel_err=" << std::abs(var_emp-var_th)/var_th   << "\n";
  return 0;
}
