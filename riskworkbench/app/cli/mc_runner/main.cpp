#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/models/gbm.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/pricing/mc_pricer.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>

static void print_usage(const char* prog) {
  std::cerr << "Usage:\n  " << prog
            << " S0 K r q sigma T seed n_target batch tol n_steps [--put] [--log] [--csv FILE]\n";
}

int main(int argc, char** argv) {
  if (argc < 12) {
    print_usage(argv[0]);
    return 1;
  }

  // Parse positionnels
  double S0, K, r, q, sigma, T;
  unsigned long long seed;
  std::size_t n_target, batch, n_steps;
  double tol;

  try {
    S0       = std::stod(argv[1]);
    K        = std::stod(argv[2]);
    r        = std::stod(argv[3]);
    q        = std::stod(argv[4]);
    sigma    = std::stod(argv[5]);
    T        = std::stod(argv[6]);
    seed     = std::stoull(argv[7]);
    n_target = static_cast<std::size_t>(std::stoull(argv[8]));
    batch    = static_cast<std::size_t>(std::stoull(argv[9]));
    tol      = std::stod(argv[10]);
    n_steps  = static_cast<std::size_t>(std::stoull(argv[11]));
  } catch (...) {
    print_usage(argv[0]);
    return 1;
  }

  // Flags optionnels
  bool is_put = false;
  bool want_log = false;
  std::string csv_file;

  for (int i = 12; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--put") {
      is_put = true;
    } else if (arg == "--log") {
      want_log = true;
    } else if (arg == "--call") {
      is_put = false;
    } else if (arg == "--csv" && i + 1 < argc) {
      csv_file = argv[++i];
    } else if (arg.rfind("--csv=", 0) == 0) {
      csv_file = arg.substr(6);
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  try {
    // Construire les objets “maison”
    rw::market::MarketData mkt(S0, r, q);
    rw::market::OptionType type = is_put ? rw::market::OptionType::Put
                                         : rw::market::OptionType::Call;
    rw::market::Instrument inst(K, T, type);

    rw::models::GbmParams gp{r, q, sigma};
    rw::models::Gbm model(gp);

    rw::config::McConfig cfg(n_target, batch, tol, n_steps, seed);
    rw::pricing::McPricer pricer(model, cfg);
    if (want_log) pricer.enable_convergence_log(true);

    // Run MC
    const auto res = pricer.price_european(mkt, inst);

    // Affichage
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);

    std::cout << "price       : " << res.price        << "\n"
              << "std_error   : " << res.std_error    << "\n"
              << "ci_low      : " << res.ci_low       << "\n"
              << "ci_high     : " << res.ci_high      << "\n"
              << "n_effective : " << res.n_effective  << "\n"
              << "elapsed_ms  : " << res.elapsed_ms   << "\n";

    // Optionnel : dump CSV du journal de convergence
    if (!csv_file.empty() && !res.convergence_log.empty()) {
      std::ofstream ofs(csv_file);
      if (!ofs) {
        std::cerr << "Cannot open CSV file: " << csv_file << "\n";
        return 2;
      }
      ofs.setf(std::ios::fixed);
      ofs << std::setprecision(10);
      ofs << "n_cum,estimate,half_width_95\n";
      for (const auto& pt : res.convergence_log) {
        ofs << pt.n_cum << ',' << pt.estimate << ',' << pt.half_width_95 << '\n';
      }
      ofs.close();
      std::cout << "convergence_log written to: " << csv_file << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 3;
  }

  return 0;
}
