#include "rw/market/surface.hpp"     // MarketSurface + read_market_surface
#include "rw/market/iv_surface.hpp"  // IvSurface
#include <iostream>
#include <string>
#include <limits>
#include <cmath>

int main(int argc, char** argv) {
  std::string path = (argc>1 ? argv[1] : "data/market_samples/sample_market.csv");

  // Charge la surface marché (S0,r,q + rows déjà filtrées/normalisées)
  rw::market::MarketSurface ms = rw::market::read_market_surface(path);

  // Construit la surface d'IV
  rw::market::IvSurface surf(ms.S0, ms.r, ms.q);
  surf.add_quotes_from_rows(ms.rows);
  surf.build();

  std::cout << "Surface built: maturities=" << surf.maturities()
            << " kept=" << surf.quotes_kept()
            << " skipped=" << surf.quotes_skipped()
            << "\n";

  // Aperçu : IV à ATM/ITM/OTM pour quelques T
  for (const double T : {0.25, 0.5, 1.0}) {
    double iv_atm = surf.iv(ms.S0, T);
    double iv_itm = surf.iv(0.9 * ms.S0, T);
    double iv_otm = surf.iv(1.1 * ms.S0, T);
    if (std::isfinite(iv_atm)) {
      std::cout << "T=" << T
                << "  iv(ATM)=" << iv_atm
                << "  iv(0.9*S0)=" << iv_itm
                << "  iv(1.1*S0)=" << iv_otm
                << "\n";
    }
  }

  return 0;
}
