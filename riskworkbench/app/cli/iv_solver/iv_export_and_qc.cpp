// app/cli/iv_solver/iv_export_and_qc.cpp
#include "rw/market/surface.hpp"     // read_market_surface()
#include "rw/market/iv_surface.hpp"  // IvSurface
#include "rw/io/market_csv.hpp"      // read_market_meta_csv (pour le nommage)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <ctime>
#include <cctype>

static inline double df(double r,double T){ return std::exp(-r*T); }
static inline double N(double x){ return 0.5 * std::erfc(-x/std::sqrt(2.0)); }

static double bs_call(double S, double K, double r, double q, double T, double sigma){
  if (T<=0.0 || sigma<=0.0) return std::max(0.0, S*df(q,T) - K*df(r,T));
  const double v = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T)/v;
  const double d2 = d1 - v;
  return S*df(q,T)*N(d1) - K*df(r,T)*N(d2);
}
static double bs_put(double S, double K, double r, double q, double T, double sigma){
  if (T<=0.0 || sigma<=0.0) return std::max(0.0, K*df(r,T) - S*df(q,T));
  const double v = sigma*std::sqrt(T);
  const double m = std::log(S/K) + (r-q)*T;
  const double d1 = (m + 0.5*sigma*sigma*T)/v;
  const double d2 = d1 - v;
  return K*df(r,T)*N(-d2) - S*df(q,T)*N(-d1);
}

static bool call_in_bounds(double S,double K,double r,double q,double T,double C){
  const double lo = std::max(0.0, S*df(q,T) - K*df(r,T));
  const double hi = S*df(q,T);
  return (C >= lo && C <= hi);
}
static bool put_in_bounds(double S,double K,double r,double q,double T,double P){
  const double lo = std::max(0.0, K*df(r,T) - S*df(q,T));
  const double hi = K*df(r,T);
  return (P >= lo && P <= hi);
}

static std::vector<double> parse_tlist(const std::string& s){
  std::vector<double> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    try { out.push_back(std::stod(tok)); } catch(...) {}
  }
  return out;
}

static void usage(const char* argv0){
  std::cerr <<
    "Usage: " << argv0 << " -f input.csv [-o grid.csv] [--nk 41] [--kmin a] [--kmax b] [--t 0.25,0.5,1]\n"
    "Options:\n"
    "  -f / --file   CSV marche (avec meta underlying,S0,r,q)\n"
    "  -o / --out    CSV de sortie (defaut auto dans data/iv_exports/<UND>/...)\n"
    "  --nk          nb de points K (def: 41)\n"
    "  --kmin        K_min (def: 0.7*S0)\n"
    "  --kmax        K_max (def: 1.3*S0)\n"
    "  --t           liste de maturites (def: 0.25,0.5,1)\n"
    "  --linearK     grille lineaire en K (defaut: log-space en x=ln(K/S0))\n"
    "  --export-w    ajoute la colonne w=sigma^2*T dans le CSV\n";
}

struct Stat { double rel; double abs; };
static void summarize(const char* title, int priced, int used, int nan, const std::vector<Stat>& v){
  std::cout << title << " (priced=" << priced << ", used=" << used << ", nan=" << nan << ")\n";
  if (v.empty()) { std::cout << "  no rows.\n"; return; }
  std::vector<double> re(v.size());
  double mean=0.0, mx=0.0;
  for (size_t i=0;i<v.size();++i){ re[i]=v[i].rel; mean += v[i].rel; if (v[i].rel>mx) mx=v[i].rel; }
  mean /= v.size(); std::sort(re.begin(), re.end());
  auto pct = [&](double p){
    double idx = p*(re.size()-1);
    size_t lo = std::floor(idx), hi = std::ceil(idx);
    if (hi==lo) return re[lo];
    double a = idx - lo; return (1.0-a)*re[lo] + a*re[hi];
  };
  std::cout << "  rel.err mean=" << mean
            << "  p50=" << pct(0.50)
            << "  p90=" << pct(0.90)
            << "  p99=" << pct(0.99)
            << "  max=" << mx << "\n";
}

// Calendar QC: w(x,T) non-decreasing in T
static void qc_calendar_arb(const rw::market::IvSurface& surf,
                            double S0,
                            const std::vector<double>& Tlist,
                            int NX = 61,
                            double xMin = std::log(0.7), double xMax = std::log(1.3)) {
  if (Tlist.size() < 2) { std::cout << "Calendar QC: need >=2 maturities.\n"; return; }
  int violations = 0, checks = 0;
  for (size_t i=0;i+1<Tlist.size();++i){
    double T1 = Tlist[i], T2 = Tlist[i+1];
    for (int j=0;j<NX;++j){
      double u = (NX==1 ? 0.0 : double(j)/(NX-1));
      double x = (1.0-u)*xMin + u*xMax;
      double K = S0 * std::exp(x);
      double s1 = surf.iv(K, T1), s2 = surf.iv(K, T2);
      if (!(std::isfinite(s1) && std::isfinite(s2))) continue;
      double w1 = s1*s1*T1, w2 = s2*s2*T2;
      ++checks;
      if (w2 + 1e-12 < w1) ++violations;
    }
  }
  std::cout << "Calendar QC: checks=" << checks << " violations=" << violations
            << "  (" << (checks? (100.0*violations/checks):0.0) << "%)\n";
}

static std::string sanitize_name(std::string s){
  std::string out; out.reserve(s.size());
  for (unsigned char c : s) {
    if (std::isalnum(c) || c=='_' || c=='-') out.push_back(char(c));
    else if (c==' ' || c=='/' || c=='\\')    out.push_back('_');
    // autres caractères => ignorés
  }
  if (out.empty()) out = "UNKNOWN";
  return out;
}

int main(int argc, char** argv){
  std::string in_path, out_path; // si vide => on génère sous data/iv_exports/<UND>/...
  int NK = 41;
  double Kmin = std::numeric_limits<double>::quiet_NaN();
  double Kmax = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> Tlist {0.25, 0.5, 1.0};
  bool logx = true;     // grille par défaut: uniforme en x = ln(K/S0)
  bool export_w = false;

  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if ((a=="-f" || a=="--file") && i+1<argc) in_path = argv[++i];
    else if ((a=="-o" || a=="--out") && i+1<argc) out_path = argv[++i];
    else if (a=="--nk" && i+1<argc) NK = std::max(2, std::atoi(argv[++i]));
    else if (a=="--kmin" && i+1<argc) Kmin = std::atof(argv[++i]);
    else if (a=="--kmax" && i+1<argc) Kmax = std::atof(argv[++i]);
    else if (a=="--t" && i+1<argc) Tlist = parse_tlist(argv[++i]);
    else if (a=="--linearK") logx = false;
    else if (a=="--export-w") export_w = true;
    else if (a=="-h" || a=="--help") { usage(argv[0]); return 0; }
  }
  if (in_path.empty()){ usage(argv[0]); return 1; }

  // 1) Lecture marche
  rw::market::MarketSurface ms;
  try { ms = rw::market::read_market_surface(in_path); }
  catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 2; }

  // underlying pour nommage (si présent)
  std::string und = "UNKNOWN";
  if (auto meta = rw::io::read_market_meta_csv(in_path)) {
    if (!meta->underlying.empty()) und = sanitize_name(meta->underlying);
  }

  // 2) Surface IV
  rw::market::IvSurface surf(ms.S0, ms.r, ms.q);
  surf.add_quotes_from_rows(ms.rows);
  surf.build();
  std::cout << "Surface built: maturities=" << surf.maturities()
            << " kept=" << surf.quotes_kept()
            << " skipped=" << surf.quotes_skipped() << "\n";

  if (!std::isfinite(Kmin)) Kmin = 0.7 * ms.S0;
  if (!std::isfinite(Kmax)) Kmax = 1.3 * ms.S0;
  if (Kmax <= Kmin) { std::cerr << "error: Kmax<=Kmin\n"; return 3; }

  // 3) Détermination du chemin de sortie par défaut dans data/
  if (out_path.empty()) {
    // data/iv_exports/<UNDERLYING>/<FILENAME>.csv
    const std::filesystem::path base = std::filesystem::path("data") / "iv_exports" / und;

    // assemble nom de fichier: ivgrid_<UND>_<logx|linK>_NKxx_Kmin-Kmax_Ttlist_YYYYmmdd-HHMMSS.csv
    std::ostringstream tlist_ss; tlist_ss.setf(std::ios::fixed); tlist_ss<<std::setprecision(4);
    for (size_t i=0;i<Tlist.size();++i){ if(i) tlist_ss<<'-'; tlist_ss<<Tlist[i]; }

    std::time_t tt = std::time(nullptr);
    std::tm tm{};
    #ifdef _WIN32
      localtime_s(&tm, &tt);
    #else
      tm = *std::localtime(&tt);
    #endif
    std::ostringstream ts; ts<< std::put_time(&tm, "%Y%m%d-%H%M%S");

    std::ostringstream fname;
    fname << "ivgrid_" << und
          << "_" << (logx ? "logx" : "linK")
          << "_NK" << NK
          << "_K" << (long)std::llround(Kmin) << "-" << (long)std::llround(Kmax)
          << "_T" << tlist_ss.str()
          << "_" << ts.str() << ".csv";

    std::filesystem::create_directories(base);
    out_path = (base / fname.str()).string();
  }

  // 4) Export grille (T,K,iv[,w],call_price,put_price)
  {
    std::filesystem::path p(out_path);
    if (p.has_parent_path()) {
      std::error_code ec;
      std::filesystem::create_directories(p.parent_path(), ec);
    }
  }

  std::ofstream out(out_path);
  if (!out) { std::cerr << "error: cannot open " << out_path << " for write\n"; return 4; }
  out << "T,K,iv";
  if (export_w) out << ",w";
  out << ",call_price,put_price\n";

  for (double T : Tlist) {
    for (int i=0;i<NK;++i) {
      double u = (NK==1 ? 0.0 : double(i)/(NK-1));
      double K = 0.0;
      if (logx) {
        double xMin = std::log(Kmin / ms.S0);
        double xMax = std::log(Kmax / ms.S0);
        double x = (1.0-u)*xMin + u*xMax;
        K = ms.S0 * std::exp(x);
      } else {
        K = (1.0-u)*Kmin + u*Kmax;
      }
      double iv = surf.iv(K, T);
      out << T << "," << K << "," << (std::isfinite(iv) ? iv : std::numeric_limits<double>::quiet_NaN());
      if (export_w) {
        double w = (std::isfinite(iv) ? (iv*iv*T) : std::numeric_limits<double>::quiet_NaN());
        out << "," << w;
      }
      double c = std::numeric_limits<double>::quiet_NaN();
      double p = std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(iv) && iv > 0.0) {
        c = bs_call(ms.S0, K, ms.r, ms.q, T, iv);
        p = bs_put (ms.S0, K, ms.r, ms.q, T, iv);
      }
      out << "," << c << "," << p << "\n";
    }
  }
  out.close();

  std::cout << "Grid exported to: " << out_path
            << "  (T in {"; for (size_t i=0;i<Tlist.size();++i){ if(i) std::cout<<","; std::cout<<Tlist[i]; } std::cout<<"}, "
            << (logx ? "log-grid in x" : "linear grid in K")
            << ", K in [" << Kmin << "," << Kmax << "], NK=" << NK << ")\n";

  // 5) QC repricing
  std::vector<Stat> errs_all, errs_arb;
  int priced=0, used_all=0, used_arb=0, nan_all=0, nan_arb=0;

  for (const auto& r : ms.rows) {
    if (!std::isfinite(r.mid_price)) continue;
    ++priced;

    const bool ok = r.is_call
      ? call_in_bounds(ms.S0, r.K, ms.r, ms.q, r.T, r.mid_price)
      : put_in_bounds (ms.S0, r.K, ms.r, ms.q, r.T, r.mid_price);

    auto push_eval = [&](std::vector<Stat>& dst, int& used, int& nan){
      double sigma = surf.iv(r.K, r.T);
      if (!std::isfinite(sigma) || sigma <= 0.0) { ++nan; return; }
      ++used;
      double theo = r.is_call ? bs_call(ms.S0, r.K, ms.r, ms.q, r.T, sigma)
                              : bs_put (ms.S0, r.K, ms.r, ms.q, r.T, sigma);
      double eabs = std::fabs(theo - r.mid_price);
      double scale = 1.0;
      if (r.is_call) scale = std::max({1.0, std::fabs(r.mid_price), ms.S0*df(ms.q, r.T)});
      else           scale = std::max({1.0, std::fabs(r.mid_price), r.K*df(ms.r, r.T)});
      double erel = eabs / scale;
      dst.push_back({erel, eabs});
    };

    // all priced
    push_eval(errs_all, used_all, nan_all);
    // only arbitrage-OK
    if (ok) push_eval(errs_arb, used_arb, nan_arb);
  }

  summarize("QC repricing (arb-OK only)", priced, used_arb, nan_arb, errs_arb);
  summarize("QC repricing (all priced)",  priced, used_all, nan_all, errs_all);

  // 6) QC calendar (croissance de w en T)
  qc_calendar_arb(surf, ms.S0, Tlist);

  return 0;
}
