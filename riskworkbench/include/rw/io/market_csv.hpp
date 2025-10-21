#pragma once
#include <string>
#include <vector>
#include <optional>
#include <limits>

namespace rw::io {

struct QuoteRow {
  double K       = std::numeric_limits<double>::quiet_NaN();
  double T       = std::numeric_limits<double>::quiet_NaN();
  bool   is_call = true;

  // prix
  double mid_price = std::numeric_limits<double>::quiet_NaN();
  double bid       = std::numeric_limits<double>::quiet_NaN();
  double ask       = std::numeric_limits<double>::quiet_NaN();

  // IV
  double iv_mid    = std::numeric_limits<double>::quiet_NaN();
  double iv_bid    = std::numeric_limits<double>::quiet_NaN();
  double iv_ask    = std::numeric_limits<double>::quiet_NaN();
};

// Lit un CSV de quotes (cf. colonnes typiques). Normalise mid si possible.
// Filtre les lignes invalides (T<=0, K<=0, spreads incohérents, etc.).
// Retourne uniquement les lignes **valides**.
// num_ignored/warnings sont optionnels pour diagnostic.
std::vector<QuoteRow>
read_market_csv(const std::string& path,
                std::size_t* num_ignored = nullptr,
                std::vector<std::string>* warnings = nullptr);

// (Optionnel) métadonnées marché si présentes dans le CSV (underlying,S0,r,q).
struct MarketMeta {
  std::string underlying;
  double S0 = std::numeric_limits<double>::quiet_NaN();
  double r  = std::numeric_limits<double>::quiet_NaN();
  double q  = std::numeric_limits<double>::quiet_NaN();
};

// Tente d’extraire underlying,S0,r,q si le CSV les contient (en-tête ou lignes méta).
std::optional<MarketMeta> read_market_meta_csv(const std::string& path);

} // namespace rw::io
