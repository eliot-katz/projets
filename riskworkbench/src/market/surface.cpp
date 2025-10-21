#include "rw/market/surface.hpp"
#include <stdexcept>
#include <cmath>

namespace rw::market {

MarketSurface read_market_surface(const std::string& path) {
  std::size_t ignored = 0;
  std::vector<std::string> warns;
  auto rows = rw::io::read_market_csv(path, &ignored, &warns);

  auto meta = rw::io::read_market_meta_csv(path);
  if (!meta.has_value() || !std::isfinite(meta->S0) || !std::isfinite(meta->r) || !std::isfinite(meta->q)) {
    throw std::runtime_error("Le CSV ne contient pas de métadonnées marché complètes (underlying,S0,r,q).");
  }

  MarketSurface ms;
  ms.S0 = meta->S0;
  ms.r  = meta->r;
  ms.q  = meta->q;
  ms.rows = std::move(rows);
  return ms;
}

} // namespace rw::market
