#include <rw/market/market_data.hpp>
#include <rw/market/instrument.hpp>
#include <rw/models/gbm.hpp>
#include <rw/config/mc_config.hpp>
#include <rw/config/greeks_config.hpp>
#include <rw/pricing/greeks.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

static void usage(const char* prog) {
  std::cerr << "Usage:\n  " << prog
            << " S0 K r q sigma T seed n_target batch tol n_steps"
            << " --greek delta|vega|rho|theta|delta_pw|gamma_pw|vega_lrm"
            << " [--put]"
            << " [--bump-rel-s0 VAL]"
            << " [--bump-abs-sigma VAL]"
            << " [--bump-abs-r VAL]"
            << " [--bump-abs-T VAL]"
            << " [--crn] [--no-crn] [--antithetic]\n";
}

int main(int argc, char** argv) {
  if (argc < 12) { usage(argv[0]); return 1; }

  // Positionnels
  double S0,K,r,q,sigma,T; unsigned long long seed;
  std::size_t n_target,batch,n_steps; double tol;
  try {
    S0 = std::stod(argv[1]); K = std::stod(argv[2]);
    r = std::stod(argv[3]);  q = std::stod(argv[4]);
    sigma = std::stod(argv[5]); T = std::stod(argv[6]);
    seed = std::stoull(argv[7]);
    n_target = static_cast<std::size_t>(std::stoull(argv[8]));
    batch    = static_cast<std::size_t>(std::stoull(argv[9]));
    tol      = std::stod(argv[10]);
    n_steps  = static_cast<std::size_t>(std::stoull(argv[11]));
  } catch (...) { usage(argv[0]); return 1; }

  // Flags
  bool is_put=false, crn=true, anti=false;
  std::string greek;
  double bump_rel_s0=0.01, bump_abs_sigma=0.01, bump_abs_r=0.0001, bump_abs_T=1.0/365.0;

  for (int i=12;i<argc;++i) {
    std::string a = argv[i];
    if (a=="--put") is_put=true;
    else if (a=="--greek" && i+1<argc) greek = argv[++i];
    else if (a=="--bump-rel-s0" && i+1<argc) bump_rel_s0 = std::stod(argv[++i]);
    else if (a=="--bump-abs-sigma" && i+1<argc) bump_abs_sigma = std::stod(argv[++i]);
    else if (a=="--bump-abs-r" && i+1<argc) bump_abs_r = std::stod(argv[++i]);
    else if (a=="--bump-abs-T" && i+1<argc) bump_abs_T = std::stod(argv[++i]);
    else if (a=="--crn") crn=true;
    else if (a=="--no-crn") crn=false;
    else if (a=="--antithetic") anti=true;
    else { std::cerr << "Unknown arg: " << a << "\n"; usage(argv[0]); return 1; }
  }
  std::transform(greek.begin(), greek.end(), greek.begin(), ::tolower);

  // ✅ Un seul test de validation, qui inclut bien vega_lrm
  if (greek!="delta" && greek!="vega" && greek!="rho" && greek!="theta"
      && greek!="delta_pw" && greek!="gamma_pw" && greek!="vega_lrm") {
    std::cerr << "Missing/invalid --greek\n"; usage(argv[0]); return 1;
  }

  try {
    rw::market::MarketData mkt(S0, r, q);
    rw::market::OptionType typ = is_put ? rw::market::OptionType::Put
                                        : rw::market::OptionType::Call;
    rw::market::Instrument inst(K, T, typ);

    rw::models::GbmParams gp{r,q,sigma};
    rw::models::Gbm model(gp);

    rw::config::McConfig mccfg(n_target, batch, tol, n_steps, seed);

    rw::config::GreeksConfig gcfg(rw::config::Greek::Delta); // init, overwritten below
    if      (greek=="delta"   ) gcfg.greek = rw::config::Greek::Delta;
    else if (greek=="vega"    ) gcfg.greek = rw::config::Greek::Vega;
    else if (greek=="rho"     ) gcfg.greek = rw::config::Greek::Rho;
    else if (greek=="theta"   ) gcfg.greek = rw::config::Greek::Theta;
    else if (greek=="gamma_pw") gcfg.greek = rw::config::Greek::Gamma; // label interne

    gcfg.bump_rel_S0    = bump_rel_s0;
    gcfg.bump_abs_sigma = bump_abs_sigma;
    gcfg.bump_abs_r     = bump_abs_r;
    gcfg.bump_abs_T     = bump_abs_T;
    gcfg.use_crn        = crn;
    gcfg.use_antithetic = anti;
    gcfg.n_steps        = n_steps;

    // PW si besoin
    if (greek=="delta_pw" || greek=="gamma_pw") gcfg.use_pathwise = true;

    // ⚠️ Vérifie que l’ordre de ton constructeur correspond bien !
    // Si ton header déclare: GreeksEstimator(Gbm, GreeksConfig, McConfig)
    // garde cette ligne; sinon ajuste.
    rw::pricing::GreeksEstimator est(model, gcfg, mccfg);

    rw::pricing::MonteCarloResult res;
    if      (greek=="delta"   ) res = est.delta_brv(mkt, inst);
    else if (greek=="vega"    ) res = est.vega_brv (mkt, inst);
    else if (greek=="rho"     ) res = est.rho_brv  (mkt, inst);
    else if (greek=="theta"   ) res = est.theta_brv(mkt, inst);
    else if (greek=="delta_pw") res = est.delta_pathwise(mkt, inst);
    else if (greek=="gamma_pw") res = est.gamma_from_delta_pw(mkt, inst);
    else if (greek=="vega_lrm") res = est.vega_lrm(mkt, inst);

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);
    std::cout << "greek      : " << greek << "\n"
              << "estimate   : " << res.price << "\n";
    std::cout << std::setprecision(9);
    std::cout << "std_error  : " << res.std_error << "\n";
    std::cout << std::setprecision(6);
    std::cout << "ci_low     : " << res.ci_low << "\n"
              << "ci_high    : " << res.ci_high << "\n"
              << "n_effective: " << res.n_effective << "\n"
              << "elapsed_ms : " << res.elapsed_ms << "\n";
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n"; return 2;
  }
  return 0;
}
