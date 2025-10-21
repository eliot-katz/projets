#include "rw/market/surface.hpp"
#include "rw/io/market_csv.hpp"
#include "rw/calibration/calib.hpp"
#include "rw/smile/smile.hpp"
#include "rw/qc/qc.hpp"
#include "rw/pricing/implied_vol.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <filesystem>
#include <algorithm>
#include <numeric>

static std::string sanitize(std::string s){
  std::string o; o.reserve(s.size());
  for (unsigned char c: s){
    if (std::isalnum(c) || c=='_' || c=='-') o.push_back((char)c);
    else if (c==' '||c=='/'||c=='\\') o.push_back('_');
  }
  return o.empty()? "UNKNOWN": o;
}

static void usage(const char* a){
  std::cerr << "Usage: " << a << " -f market.csv [--tick 0.01] [--thresh 0.003]\n"
               "                 [--price-from-iv] [--only-consistent TOL]\n"
               "  --tick            : test stabilité IV à ±tick en prix (défaut: 0 = off)\n"
               "  --thresh          : succès si RMSE_price_abs < thresh * S0 (défaut 0.003 = 0.3%)\n"
               "  --price-from-iv   : utiliser BS(IV_mkt) comme price_mkt quand IV dispo\n"
               "  --only-consistent : filtrer |mid - BS(IV_mkt)| ≤ TOL (en prix)\n";
}

// BS prix C/P (local)
static inline double N(double x){ return 0.5 * std::erfc(-x * M_SQRT1_2); }
static double bs_price_cp(double S0,double K,double r,double q,double T,double sigma,bool is_call){
  if (!(std::isfinite(S0)&&std::isfinite(K)&&std::isfinite(r)&&std::isfinite(q)&&std::isfinite(T)&&std::isfinite(sigma))
      || S0<=0.0 || K<=0.0 || T<0.0) return std::numeric_limits<double>::quiet_NaN();
  if (T==0.0 || sigma<=0.0){
    const double df_r=std::exp(-r*T), df_q=std::exp(-q*T);
    const double fwd=S0*df_q - K*df_r;
    return is_call? std::max(0.0,fwd) : std::max(0.0,-fwd);
  }
  const double df_r=std::exp(-r*T), df_q=std::exp(-q*T);
  const double vs = sigma*std::sqrt(T);
  const double m  = std::log(S0/K) + (r - q + 0.5*sigma*sigma)*T;
  const double d1 = m / vs;
  const double d2 = d1 - vs;
  return is_call ? S0*df_q*N(d1) - K*df_r*N(d2)
                 : K*df_r*N(-d2) - S0*df_q*N(-d1);
}

int main(int argc, char** argv){
  std::string path;
  double tick = 0.0;
  double thresh = 0.003;      // 0.3% S0
  bool price_from_iv = false;
  double only_consistent_tol = -1.0; // ≤0 => off

  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if ((a=="-f"||a=="--file") && i+1<argc) path = argv[++i];
    else if (a=="--tick" && i+1<argc) tick = std::atof(argv[++i]);
    else if (a=="--thresh" && i+1<argc) thresh = std::atof(argv[++i]);
    else if (a=="--price-from-iv") price_from_iv = true;
    else if (a=="--only-consistent" && i+1<argc) only_consistent_tol = std::atof(argv[++i]);
    else if (a=="-h"||a=="--help"){ usage(argv[0]); return 0; }
  }
  if (path.empty()) { usage(argv[0]); return 1; }

  rw::market::MarketSurface ms;
  try { ms = rw::market::read_market_surface(path); }
  catch (const std::exception& e){ std::cerr<<"error: "<<e.what()<<"\n"; return 2; }

  std::string und = "UNKNOWN";
  if (auto meta = rw::io::read_market_meta_csv(path)) if (!meta->underlying.empty()) und = sanitize(meta->underlying);

  auto smile = rw::calib::build_smile_surface(ms);
  auto rows  = rw::qc::compare_fit(ms, smile);
  if (rows.empty()){ std::cerr<<"No comparable rows.\n"; return 0; }

  struct RowX {
    rw::qc::ErrorRow r;
    double price_mkt_iv{std::numeric_limits<double>::quiet_NaN()};
    double price_consistency_diff{std::numeric_limits<double>::quiet_NaN()}; // mid - BS(IV_mkt)
    bool keep{true};
  };
  std::vector<RowX> rx; rx.reserve(rows.size());

  for (auto r : rows){
    // prix cohérent IV marché (utilise le VRAI type call/put)
    double price_iv = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(r.iv_mkt) && r.iv_mkt>0.0)
      price_iv = bs_price_cp(ms.S0, r.K, ms.r, ms.q, r.T, r.iv_mkt, r.is_call);

    double diff = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(price_iv) && std::isfinite(r.price_mkt))
      diff = r.price_mkt - price_iv;

    bool keep = true;
    if (only_consistent_tol > 0.0 && std::isfinite(diff))
      keep = (std::fabs(diff) <= only_consistent_tol);

    if (price_from_iv && std::isfinite(price_iv)){
      r.price_mkt = price_iv;
      r.err_price = r.price_fit - r.price_mkt;
    }

    rx.push_back(RowX{r, price_iv, diff, keep});
  }

  auto compute_stats = [&](const std::vector<RowX>& v, bool only_keep){
    double sse_p=0.0, sse_iv=0.0; std::size_t n=0;
    std::vector<double> errs; errs.reserve(v.size());
    for (const auto& e : v){
      if (only_keep && !e.keep) continue;
      sse_p  += e.r.err_price*e.r.err_price;
      sse_iv += e.r.err_iv   *e.r.err_iv;
      errs.push_back(e.r.err_price);
      ++n;
    }
    if (n==0) return std::tuple<double,double,double,std::size_t,std::size_t>(NAN,NAN,NAN,0,0);
    const double rmse_price = std::sqrt(sse_p / n);
    const double rmse_iv    = std::sqrt(sse_iv / n);
    double mean_e=0.0; for (double z: errs) mean_e += z; mean_e/=n;
    double var_e=0.0; for (double z: errs) var_e += (z-mean_e)*(z-mean_e); var_e = (n>1? var_e/(n-1):0.0);
    const double sd_e = std::sqrt(var_e);
    std::size_t out3 = 0; if (sd_e>0) for (double z: errs) if (std::fabs((z-mean_e)/sd_e)>3.0) ++out3;
    return std::tuple<double,double,double,std::size_t,std::size_t>(rmse_price, rmse_iv, sd_e, n, out3);
  };

  auto [rmse_all, rmse_iv_all, sd_all, n_all, out3_all] = compute_stats(rx, false);
  auto [rmse_con, rmse_iv_con, sd_con, n_con, out3_con] = compute_stats(rx, (only_consistent_tol>0.0));

  const bool pass_all = std::isfinite(rmse_all) && (rmse_all/ms.S0 <= thresh);
  const bool pass_con = (only_consistent_tol>0.0 && std::isfinite(rmse_con)) ? (rmse_con/ms.S0 <= thresh) : false;

  // Export
  std::filesystem::path outdir = std::filesystem::path("data") / "iv_exports" / und;
  std::filesystem::create_directories(outdir);
  std::filesystem::path outcsv = outdir / "qc_residuals.csv";
  std::ofstream f(outcsv);
  f << "K,T,is_call,price_mkt,price_fit,err_price,iv_mkt,iv_fit,err_iv,price_mkt_iv,price_consistency_diff,keep\n";
  for (const auto& e: rx){
    f << e.r.K << "," << e.r.T << ","
      << (e.r.is_call?1:0) << ","
      << e.r.price_mkt << "," << e.r.price_fit << "," << e.r.err_price << ","
      << e.r.iv_mkt << "," << e.r.iv_fit << "," << e.r.err_iv << ","
      << e.price_mkt_iv << "," << e.price_consistency_diff << ","
      << (e.keep?1:0) << "\n";
  }
  f.close();

  // Stabilité (optionnelle)
  double sens_med = std::numeric_limits<double>::quiet_NaN();
  if (tick>0.0){
    std::vector<double> sens; sens.reserve(rx.size());
    for (const auto& e: rx){
      const auto& r = e.r;
      auto up = rw::bs::implied_vol_cp(ms.S0, r.K, ms.r, ms.q, r.T, r.price_mkt + tick, r.is_call);
      auto dn = rw::bs::implied_vol_cp(ms.S0, r.K, ms.r, ms.q, r.T, r.price_mkt - tick, r.is_call);
      if (up.converged && dn.converged && std::isfinite(up.sigma) && std::isfinite(dn.sigma))
        sens.push_back( (up.sigma - dn.sigma) / (2.0*tick) );
    }
    if (!sens.empty()){
      std::sort(sens.begin(), sens.end());
      sens_med = sens[sens.size()/2];
    }
  }

  std::cout.setf(std::ios::fixed); std::cout.precision(8);
  std::cout << "QC (all rows): n=" << n_all
            << "  RMSE_price=" << rmse_all
            << "  RMSE_price_rel_S0=" << (rmse_all/ms.S0)*100.0 << "%  (thresh="<< (thresh*100.0) <<"%)  "
            << (pass_all? "PASS":"FAIL") << "\n"
            << "    RMSE_iv=" << rmse_iv_all
            << "    outliers |z|>3: " << out3_all << "\n";

  if (only_consistent_tol>0.0) {
    std::cout << "QC (only consistent: |mid - BS(IV_mkt)| ≤ " << only_consistent_tol << "): n=" << n_con
              << "  RMSE_price=" << rmse_con
              << "  RMSE_price_rel_S0=" << (rmse_con/ms.S0)*100.0 << "%  "
              << (pass_con? "PASS":"FAIL") << "\n";
  }

  if (tick>0.0)
    std::cout << "Stability: median dSigma/dPrice @tick=" << tick << " = " << sens_med << "\n";

  std::cout << "CSV: \"" << outcsv.string() << "\"\n";
  return 0;
}
