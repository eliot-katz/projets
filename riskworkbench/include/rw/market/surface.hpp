#pragma once
#include <vector>
#include <string>
#include "rw/io/market_csv.hpp"

namespace rw::market {

struct MarketSurface {
  double S0 = 0.0;
  double r  = 0.0;
  double q  = 0.0;
  std::vector<rw::io::QuoteRow> rows;
};

// Construit une MarketSurface depuis un CSV.
// - lit les quotes (normalisées+filtrées)
// - lit S0,r,q si disponibles (sinon lève si absents)
MarketSurface read_market_surface(const std::string& path);

} // namespace rw::market
