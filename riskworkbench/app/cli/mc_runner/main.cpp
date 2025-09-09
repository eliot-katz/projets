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
#include <algorithm>

static void print_usage(const char* prog) {
  std::cerr << "Usage:\n  " << prog
            << " S0 K r q sigma T seed n_target batch tol n_steps"
            << " [--put] [--log] [--antithetic]"
            << " [--cv] [--cv-custom-K VAL] [--cv-type call|put]"
            << " [--cv-beta-fixed VAL] [--cv-pilot N] [--csv FILE]\n";
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
  bool use_antithetic = false;

  bool use_cv = false;
  rw::config::CvMode cv_mode = rw::config::CvMode::SameInstrument;
  double cvK_override = -1.0; // <0 means not set
  rw::market::OptionType cv_type = rw::market::OptionType::Call; // placeholder, will be synced later
  bool cv_type_set = false;

  rw::config::BetaMode beta_mode = rw::config::BetaMode::Pilot;
  double beta_fixed = 0.0;
  std::size_t cv_pilot_paths = 20000;

  std::string csv_file;

  for (int i = 12; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--put") {
      is_put = true;
    } else if (arg == "--log") {
      want_log = true;
    } else if (arg == "--antithetic") {
      use_antithetic = true;
    } else if (arg == "--cv") {
      use_cv = true;
      cv_mode = rw::config::CvMode::SameInstrument;
    } else if (arg == "--cv-custom-K" && i + 1 < argc) {
      use_cv = true;
      cv_mode = rw::config::CvMode::CustomInstrument;
      cvK_override = std::stod(argv[++i]);
    } else if (arg == "--cv-type" && i + 1 < argc) {
      std::string t = argv[++i];
      std::transform(t.begin(), t.end(), t.begin(), ::tolower);
      if (t == "call") cv_type = rw::market::OptionType::Call;
      else if (t == "put") cv_type = rw::market::OptionType::Put;
      else { std::cerr << "Unknown --cv-type: " << t << "\n"; return 1; }
      cv_type_set = true;
      use_cv = true;
    } else if (arg == "--cv-beta-fixed" && i + 1 < argc) {
      beta_mode = rw::config::BetaMode::Fixed;
      beta_fixed = std::stod(argv[++i]);
      use_cv = true;
    } else if (arg == "--cv-pilot" && i + 1 < argc) {
      cv_pilot_paths = static_cast<std::size_t>(std::stoull(argv[++i]));
      use_cv = true;
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
    // Données de marché et instrument cible
    rw::market::MarketData mkt(S0, r, q);
    rw::market::OptionType type = is_put ? rw::market::OptionType::Put
                                         : rw::market::OptionType::Call;
    rw::market::Instrument inst(K, T, type);

    // Modèle & config
    rw::models::GbmParams gp{r, q, sigma};
    rw::models::Gbm model(gp);

    // Construire un cv_instrument par défaut (sera ignoré si SameInstrument)
    rw::market::OptionType cv_eff_type = cv_type_set ? cv_type : type; // par défaut même type que la cible
    double cvK = (cvK_override > 0.0) ? cvK_override : K; // par défaut même strike
    rw::market::Instrument cv_ins(cvK, T, cv_eff_type);

    rw::config::McConfig cfg(
      n_target, batch, tol, n_steps, seed,
      use_antithetic,
      use_cv,
      cv_mode,
      cv_ins,
      beta_mode,
      beta_fixed,
      cv_pilot_paths
    );

    rw::pricing::McPricer pricer(model, cfg);
    if (want_log) pricer.enable_convergence_log(true);

    // Run MC
    const auto res = pricer.price_european(mkt, inst);

    // Affichage
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);

    const double half_width = 0.5 * (res.ci_high - res.ci_low);

    std::cout << "price         : " << res.price        << "\n"
              << "std_error     : " << res.std_error    << "\n"
              << "ci_low        : " << res.ci_low       << "\n"
              << "ci_high       : " << res.ci_high      << "\n"
              << "half_width    : " << half_width       << "\n"
              << "n_effective   : " << res.n_effective  << "\n"
              << "elapsed_ms    : " << res.elapsed_ms   << "\n"
              << "antithetic    : " << (use_antithetic ? "true" : "false") << "\n"
              << "control_var   : " << (use_cv ? "true" : "false") << "\n";

    if (use_cv) {
      std::cout << "cv_mode       : " << (cv_mode == rw::config::CvMode::SameInstrument ? "SameInstrument" : "CustomInstrument") << "\n"
                << "cv_beta_mode  : " << (beta_mode == rw::config::BetaMode::Pilot ? "Pilot" : "Fixed") << "\n";
      if (beta_mode == rw::config::BetaMode::Fixed) {
        std::cout << "cv_beta_fixed : " << beta_fixed << "\n";
      }
      if (cv_mode == rw::config::CvMode::CustomInstrument) {
        std::cout << "cv_instrument : "
                  << (cv_eff_type == rw::market::OptionType::Call ? "Call" : "Put")
                  << "(K=" << cvK << ", T=" << T << ")\n";
      }
    }

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
      std::cout << "convergence_log written to: " << csv_file
                << " (n_cum counts *physical paths*, includes only main phase)\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 3;
  }

  return 0;
}
