#include "rw/market/surface.hpp"
#include "rw/io/market_csv.hpp"
#include "rw/pricing/implied_vol.hpp"
#include "rw/calibration/svi.hpp"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

static std::string sanitize(std::string s){
  std::string o; o.reserve(s.size());
  for (unsigned char c: s){
    if (std::isalnum(c) || c=='_' || c=='-') o.push_back((char)c);
    else if (c==' '||c=='/'||c=='\\') o.push_back('_');
  }
  return o.empty()? "UNKNOWN": o;
}

int main(int argc, char** argv){
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " -f market.csv [--tolT 1e-6]\n";
    return 1;
  }
  std::string path; double tolT = 1e-6;
  for (int i=1;i<argc;++i){
    std::string a = argv[i];
    if ((a=="-f"||a=="--file") && i+1<argc) path = argv[++i];
    else if (a=="--tolT" && i+1<argc) tolT = std::atof(argv[++i]);
    else if (a=="-h"||a=="--help"){ std::cerr<<"Usage: "<<argv[0]<<" -f market.csv [--tolT 1e-6]\n"; return 0; }
  }
  if (path.empty()) { std::cerr<<"missing -f\n"; return 1; }

  rw::market::MarketSurface ms;
  try { ms = rw::market::read_market_surface(path); }
  catch (const std::exception& e){ std::cerr<<"error: "<<e.what()<<"\n"; return 2; }

  std::string und = "UNKNOWN";
  if (auto meta = rw::io::read_market_meta_csv(path)) if (!meta->underlying.empty()) und = sanitize(meta->underlying);

  // Regroupe par maturité (tolérance |ΔT|<tolT)
  std::map<double, std::vector<const rw::io::QuoteRow*>> groups;
  auto keyT = [&](double T)->double{
    // quantize by tolT
    if (!(std::isfinite(T))) return std::numeric_limits<double>::quiet_NaN();
    long long k = std::llround(T / std::max(tolT, 1e-12));
    return k * std::max(tolT, 1e-12);
  };
  for (const auto& r : ms.rows) {
    if (!(std::isfinite(r.K) && std::isfinite(r.T)) || r.K<=0.0 || r.T<=0.0) continue;
    groups[keyT(r.T)].push_back(&r);
  }

  // dossiers export
  std::filesystem::path base = std::filesystem::path("data") / "iv_exports" / und;
  std::filesystem::create_directories(base);

  std::ofstream pfile(base / "svi_params.csv");
  pfile << "T,a,b,rho,m,sigma,rmse_w,n,converged\n";

  rw::svi::Surface sviSurf;
  sviSurf.S0 = ms.S0; sviSurf.r = ms.r; sviSurf.q = ms.q;

  int total_slices=0, ok_slices=0;
  for (auto& [Tk, vec] : groups) {
    if (vec.size() < 3) continue;
    ++total_slices;

    // construit {k, w}
    std::vector<double> kx, wv;
    kx.reserve(vec.size()); wv.reserve(vec.size());
    for (auto* pr : vec) {
      double iv = std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(pr->iv_mid) && pr->iv_mid>0.0) iv = pr->iv_mid;
      else if (std::isfinite(pr->mid_price)) {
        auto res = rw::bs::implied_vol_cp(ms.S0, pr->K, ms.r, ms.q, pr->T, pr->mid_price, pr->is_call);
        if (res.converged) iv = res.sigma;
      }
      if (!std::isfinite(iv) || iv<=0.0) continue;
      double k = std::log(pr->K / ms.S0);
      double w = iv*iv * pr->T;
      kx.push_back(k); wv.push_back(w);
    }
    if (kx.size() < 3) continue;

    // ordonne par k
    std::vector<std::size_t> idx(kx.size());
    iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto i, auto j){ return kx[i] < kx[j]; });

    std::vector<double> k_sorted, w_sorted;
    k_sorted.reserve(idx.size()); w_sorted.reserve(idx.size());
    for (auto i : idx){ k_sorted.push_back(kx[i]); w_sorted.push_back(wv[i]); }

    rw::svi::Params p;
    auto rep = rw::svi::fit_slice(k_sorted, w_sorted, p);

    pfile << Tk << "," << p.a << "," << p.b << "," << p.rho << "," << p.m << "," << p.sigma
          << "," << rep.rmse_w << "," << rep.n << "," << (rep.converged?1:0) << "\n";

    rw::svi::Slice sl; sl.T = Tk; sl.p = p;
    sviSurf.slices.push_back(sl);

    if (rep.converged) ++ok_slices;

    // export smile reconstitué sur une grille de K
    const int NK = 61;
    double Kmin = 0.7*ms.S0, Kmax = 1.3*ms.S0;
    std::ofstream sf(base / ("svi_smile_T=" + std::to_string(Tk) + ".csv"));
    sf << "K,iv_fit,w_fit\n";
    for (int i=0;i<NK;++i){
      double u = (double)i/(NK-1);
      double K = std::exp( (1.0-u)*std::log(Kmin) + u*std::log(Kmax) ); // log-grid
      double k = std::log(K/ms.S0);
      double w = rw::svi::w_svi(k, p);
      double iv = (std::isfinite(w) && w>0.0) ? std::sqrt(w / Tk) : std::numeric_limits<double>::quiet_NaN();
      sf << K << "," << iv << "," << w << "\n";
    }
  }

  std::sort(sviSurf.slices.begin(), sviSurf.slices.end(),
            [](const rw::svi::Slice& a, const rw::svi::Slice& b){ return a.T < b.T; });

  std::cout << "SVI: slices=" << total_slices << " converged=" << ok_slices << "\n";
  std::cout << "Params file: " << (base / "svi_params.csv") << "\n";
  std::cout << "Smiles CSV:  " << (base) << "/svi_smile_T=*.csv\n";

  // petit QC ponctuel
  if (!sviSurf.slices.empty()){
    double Tmid = sviSurf.slices[sviSurf.slices.size()/2].T;
    double K_atm = ms.S0;
    double iv_atm = sviSurf.iv(K_atm, Tmid);
    std::cout << "QC: iv(K=S0,T~mid) = " << iv_atm << "\n";
  }

  return 0;
}
