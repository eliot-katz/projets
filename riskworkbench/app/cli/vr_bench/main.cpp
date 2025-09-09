#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/models/gbm.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/pricing/mc_pricer.hpp>
#include <rw/pricing/analytic_bs.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

struct Case {
  double S0, K, r, q, sigma, T;
  rw::market::OptionType type;
  const char* name;
};

static void print_usage(const char* prog) {
  std::cerr << "Usage:\n  " << prog
            << " S0 K r q sigma T seed N batch n_steps [--put] [--pilot Npilot]\n"
            << "If no args are given, runs two defaults: ATM 1Y (call) and OTM 3M (call).\n";
}

struct RunOut {
  double price;
  double se;
  double lo;
  double hi;
  std::size_t n_eff;
};

static RunOut run_mode(const Case& c,
                       unsigned long long seed,
                       std::size_t N,
                       std::size_t batch,
                       std::size_t n_steps,
                       bool use_antithetic,
                       bool use_cv,
                       std::size_t pilot_paths)
{
  rw::market::MarketData mkt(c.S0, c.r, c.q);
  rw::market::Instrument inst(c.K, c.T, c.type);
  rw::models::Gbm model({c.r, c.q, c.sigma});

  // CV config: SameInstrument + Pilot
  rw::config::McConfig cfg(
    N, batch, -1.0, n_steps, seed,
    use_antithetic,
    use_cv,
    rw::config::CvMode::SameInstrument,
    rw::market::Instrument(c.K, c.T, c.type), // ignored when SameInstrument
    rw::config::BetaMode::Pilot,
    0.0,
    pilot_paths
  );

  rw::pricing::McPricer pricer(model, cfg);
  const auto res = pricer.price_european(mkt, inst);
  return {res.price, res.std_error, res.ci_low, res.ci_high, res.n_effective};
}

static double bs_ref_price(const Case& c) {
  if (c.type == rw::market::OptionType::Call) {
    return rw::pricing::price_call_bs(c.S0, c.K, c.r, c.q, c.sigma, c.T);
  } else {
    return rw::pricing::price_put_bs (c.S0, c.K, c.r, c.q, c.sigma, c.T);
  }
}

static void bench_one(const Case& c,
                      unsigned long long seed,
                      std::size_t N,
                      std::size_t batch,
                      std::size_t n_steps,
                      std::size_t pilot_paths)
{
  const double bs = bs_ref_price(c);

  const RunOut plain   = run_mode(c, seed, N, batch, n_steps, /*anti*/false, /*cv*/false, pilot_paths);
  const RunOut anti    = run_mode(c, seed, N, batch, n_steps, /*anti*/true,  /*cv*/false, pilot_paths);
  const RunOut cv      = run_mode(c, seed, N, batch, n_steps, /*anti*/false, /*cv*/true,  pilot_paths);
  const RunOut anti_cv = run_mode(c, seed, N, batch, n_steps, /*anti*/true,  /*cv*/true,  pilot_paths);

  auto in_ci = [&](const RunOut& r){ return (bs >= r.lo && bs <= r.hi); };

  const double R_anti    = (plain.se > 0 && anti.se > 0)    ? std::pow(plain.se / anti.se, 2.0)    : 0.0;
  const double R_cv      = (plain.se > 0 && cv.se > 0)      ? std::pow(plain.se / cv.se, 2.0)      : 0.0;
  const double R_anti_cv = (plain.se > 0 && anti_cv.se > 0) ? std::pow(plain.se / anti_cv.se, 2.0) : 0.0;

  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(6);

  std::cout << "\n=== VR Bench: " << c.name << " ("
            << (c.type == rw::market::OptionType::Call ? "Call" : "Put")
            << ", S0="<<c.S0<<", K="<<c.K<<", r="<<c.r<<", q="<<c.q
            << ", sigma="<<c.sigma<<", T="<<c.T<<") ===\n";
  std::cout << "seed="<<seed<<", N="<<N<<" (physical paths), batch="<<batch<<", steps="<<n_steps
            << ", pilot_paths="<<pilot_paths<<" (used only for CV)\n";
  std::cout << "BS_ref=" << bs << "\n\n";

  std::cout << std::left
            << std::setw(12) << "Mode"
            << std::setw(14) << "price"
            << std::setw(14) << "std_error"
            << std::setw(14) << "ci_low"
            << std::setw(14) << "ci_high"
            << std::setw(12) << "in_BS_CI"
            << std::setw(10) << "R(wrt P)"
            << "\n";
  std::cout << std::string(90, '-') << "\n";

  auto row = [&](const char* lbl, const RunOut& r, double R){
    std::cout << std::left
              << std::setw(12) << lbl
              << std::setw(14) << r.price
              << std::setw(14) << r.se
              << std::setw(14) << r.lo
              << std::setw(14) << r.hi
              << std::setw(12) << (in_ci(r) ? "YES" : "NO")
              << std::setw(10) << (R > 0 ? R : 1.0)
              << "\n";
  };

  row("plain",    plain,    1.0);
  row("anti",     anti,     R_anti);
  row("cv",       cv,       R_cv);
  row("anti+cv",  anti_cv,  R_anti_cv);

  // Small summary / checks
  const bool var_ok = (R_anti >= 1.0 - 1e-9) && (R_cv >= 1.0 - 1e-9) && (R_anti_cv >= 1.0 - 1e-9);
  const bool cv_better_or_eq = (anti_cv.se <= cv.se + 1e-12);

  std::cout << "\nChecks: ";
  std::cout << "[prices ok? all in CI] "
            << ((in_ci(plain) && in_ci(anti) && in_ci(cv) && in_ci(anti_cv)) ? "PASS" : "FAIL");
  std::cout << " | [variance ↓? R>=1] " << (var_ok ? "PASS" : "FAIL");
  std::cout << " | [anti+cv <= cv (SE)] " << (cv_better_or_eq ? "PASS" : "FAIL");
  std::cout << "\n";
}

int main(int argc, char** argv) {
  if (argc == 1) {
    // Defaults: ATM 1Y (call) and OTM 3M (call)
    const unsigned long long seed = 42ULL;
    const std::size_t N = 400000;
    const std::size_t batch = 100000;
    const std::size_t steps = 1;
    const std::size_t pilot_paths = 20000;

    std::vector<Case> cases = {
      {100.0, 100.0, 0.02, 0.00, 0.20, 1.00, rw::market::OptionType::Call, "ATM 1Y"},
      {100.0, 110.0, 0.02, 0.00, 0.20, 0.25, rw::market::OptionType::Call, "OTM 3M"}
    };

    for (const auto& c : cases) bench_one(c, seed, N, batch, steps, pilot_paths);
    return 0;
  }

  if (argc < 12) {
    print_usage(argv[0]);
    return 1;
  }

  // Parse a single case from CLI
  double S0, K, r, q, sigma, T;
  unsigned long long seed;
  std::size_t N, batch, steps;
  std::size_t pilot_paths = 20000;
  rw::market::OptionType type = rw::market::OptionType::Call;

  try {
    S0    = std::stod(argv[1]);
    K     = std::stod(argv[2]);
    r     = std::stod(argv[3]);
    q     = std::stod(argv[4]);
    sigma = std::stod(argv[5]);
    T     = std::stod(argv[6]);
    seed  = std::stoull(argv[7]);
    N     = static_cast<std::size_t>(std::stoull(argv[8]));
    batch = static_cast<std::size_t>(std::stoull(argv[9]));
    steps = static_cast<std::size_t>(std::stoull(argv[10]));
  } catch (...) {
    print_usage(argv[0]);
    return 1;
  }

  for (int i = 11; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--put") {
      type = rw::market::OptionType::Put;
    } else if (a == "--pilot" && i + 1 < argc) {
      pilot_paths = static_cast<std::size_t>(std::stoull(argv[++i]));
    } else {
      std::cerr << "Unknown option: " << a << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  Case c{S0, K, r, q, sigma, T, type, "custom"};
  bench_one(c, seed, N, batch, steps, pilot_paths);
  return 0;
}
