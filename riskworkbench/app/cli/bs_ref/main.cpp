#include <rw/pricing/analytic_bs.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>

struct Case {
  double S0, K, r, q, sigma, T;
};

static void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " [S0 K r q sigma T]\n"
            << "If no arguments are provided, runs 3 reference cases.\n";
}

int main(int argc, char** argv) {
  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(6);

  std::vector<Case> cases;
  if (argc == 1) {
    // 2–3 cas “or” pour valider rapidement
    cases.push_back({100.0, 100.0, 0.02, 0.00, 0.20, 1.00});
    cases.push_back({100.0, 110.0, 0.02, 0.00, 0.25, 0.50});
    cases.push_back({80.0,  100.0, -0.01, 0.01, 0.35, 2.00}); // r/q potentiellement négatifs
  } else if (argc == 7) {
    Case c;
    try {
      c.S0    = std::stod(argv[1]);
      c.K     = std::stod(argv[2]);
      c.r     = std::stod(argv[3]);
      c.q     = std::stod(argv[4]);
      c.sigma = std::stod(argv[5]);
      c.T     = std::stod(argv[6]);
    } catch (...) {
      print_usage(argv[0]);
      return 1;
    }
    cases.push_back(c);
  } else {
    print_usage(argv[0]);
    return 1;
  }

  std::cout << " S0       K        r        q     sigma        T        Call        Put    ParityGap\n";
  std::cout << "-------------------------------------------------------------------------------------\n";
  for (const auto& c : cases) {
    const double call = rw::pricing::price_call_bs(c.S0, c.K, c.r, c.q, c.sigma, c.T);
    const double put  = rw::pricing::price_put_bs (c.S0, c.K, c.r, c.q, c.sigma, c.T);
    const double gap  = rw::pricing::put_call_parity_gap(call, put, c.S0, c.K, c.r, c.q, c.T);

    std::cout << std::setw(7)  << c.S0   << ' '
              << std::setw(7)  << c.K    << ' '
              << std::setw(8)  << c.r    << ' '
              << std::setw(8)  << c.q    << ' '
              << std::setw(8)  << c.sigma<< ' '
              << std::setw(9)  << c.T    << ' '
              << std::setw(11) << call   << ' '
              << std::setw(10) << put    << ' '
              << std::setw(11) << gap    << '\n';
  }
  return 0;
}
