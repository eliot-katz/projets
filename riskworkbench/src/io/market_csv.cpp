#include "rw/io/market_csv.hpp"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace {

// --- helpers texte ---
static inline std::string trim(std::string s) {
  auto notsp = [](int ch){ return !std::isspace(ch); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), notsp));
  s.erase(std::find_if(s.rbegin(), s.rend(), notsp).base(), s.end());
  return s;
}
static inline std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
  return s;
}

// CSV splitter minimal qui gère les champs entre "..."
static std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::string field;
  bool in_quotes = false;
  for (size_t i=0;i<line.size();++i) {
    char c = line[i];
    if (c == '"') {
      if (in_quotes && i+1<line.size() && line[i+1] == '"') { field.push_back('"'); ++i; }
      else { in_quotes = !in_quotes; }
    } else if (c == ',' && !in_quotes) {
      out.push_back(trim(field)); field.clear();
    } else {
      field.push_back(c);
    }
  }
  out.push_back(trim(field));
  return out;
}

static inline bool is_nan(double x){ return std::isnan(x); }
static inline bool is_finite(double x){ return std::isfinite(x); }

// parse double tolérant (“” -> NaN)
static double parse_double(const std::string& s) {
  if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
  char* end=nullptr;
  double v = std::strtod(s.c_str(), &end);
  if (end==s.c_str()) return std::numeric_limits<double>::quiet_NaN();
  return v;
}

static bool parse_is_call(std::string s) {
  s = lower(trim(s));
  if (s=="c" || s=="call" || s=="1" || s=="true" || s=="t") return true;
  if (s=="p" || s=="put"  || s=="0" || s=="false"|| s=="f") return false;
  // défaut: call
  return true;
}

// récupère index de colonne via map (synonymes acceptés)
static int col(const std::unordered_map<std::string,int>& idx, std::initializer_list<const char*> names) {
  for (auto* n: names) {
    auto it = idx.find(lower(n));
    if (it != idx.end()) return it->second;
  }
  return -1;
}

static bool invalid_price_spread(double bid, double ask) {
  if (is_nan(bid) || is_nan(ask)) return false; // insuffisant pour juger
  if (ask < bid) return true;                   // incohérent
  double mid = 0.5*(bid+ask);
  if (mid <= 0.0) return false;
  double rel = (ask-bid)/mid;
  return rel > 0.6; // “aberrant”: spread > 60% du mid (ajuste au besoin)
}

static bool invalid_iv_spread(double iv_bid, double iv_ask) {
  if (is_nan(iv_bid) || is_nan(iv_ask)) return false;
  if (iv_ask < iv_bid) return true;
  if (iv_bid <= 0.0 || iv_ask <= 0.0) return true;
  double mid = 0.5*(iv_bid+iv_ask);
  double rel = (iv_ask-iv_bid)/mid;
  return (iv_bid <= 0.0) || (iv_ask <= 0.0) || rel > 1.0 || iv_ask > 5.0; // IV > 500% jugée aberrante
}

} // namespace

namespace rw::io {

std::optional<MarketMeta> read_market_meta_csv(const std::string& path) {
  std::ifstream f(path);
  if (!f) return std::nullopt;

  std::string line;
  std::vector<std::string> header;
  bool header_seen = false;
  std::optional<MarketMeta> meta;
  while (std::getline(f, line)) {
    if (!line.empty() && line.back()=='\r') line.pop_back();
    auto l = trim(line);
    if (l.empty() || l.rfind("#",0)==0) continue;

    auto cells = split_csv_line(l);

    if (!header_seen) {
      header = cells;
      header_seen = true;
      for (auto& h: header) h = lower(trim(h));
      continue;
    }

    // map colonnes
    std::unordered_map<std::string,int> idx;
    for (int i=0;i<(int)header.size();++i) idx[header[i]] = i;

    int iS0 = col(idx, {"s0"});
    int ir  = col(idx, {"r","rate"});
    int iq  = col(idx, {"q","div","dividend"});
    int iUnd= col(idx, {"underlying","asset"});

    // S’il manque la grille K/T on peut être sur une ligne “métadonnée” (ou la 1ère avec tout)
    if (iS0>=0 || ir>=0 || iq>=0 || iUnd>=0) {
      MarketMeta m;
      if (iUnd>=0 && iUnd < (int)cells.size()) m.underlying = cells[iUnd];
      if (iS0>=0  && iS0  < (int)cells.size()) m.S0 = parse_double(cells[iS0]);
      if (ir>=0   && ir   < (int)cells.size()) m.r  = parse_double(cells[ir]);
      if (iq>=0   && iq   < (int)cells.size()) m.q  = parse_double(cells[iq]);
      meta = m;
      break; // on prend la première occurrence
    }
  }
  return meta;
}

std::vector<QuoteRow>
read_market_csv(const std::string& path,
                std::size_t* num_ignored,
                std::vector<std::string>* warnings)
{
  if (num_ignored) *num_ignored = 0;
  std::vector<QuoteRow> out;

  std::ifstream f(path);
  if (!f) {
    if (warnings) warnings->push_back("Impossible d'ouvrir le fichier: " + path);
    return out;
  }

  std::string line;
  std::vector<std::string> header;
  bool header_seen = false;

  while (std::getline(f, line)) {
    if (!line.empty() && line.back()=='\r') line.pop_back();
    auto l = trim(line);
    if (l.empty() || l.rfind("#",0)==0) continue;

    auto cells = split_csv_line(l);

    if (!header_seen) {
      header = cells;
      header_seen = true;
      for (auto& h: header) h = lower(trim(h));
      continue;
    }

    // index colonnes (synonymes autorisés)
    std::unordered_map<std::string,int> idx;
    for (int i=0;i<(int)header.size();++i) idx[header[i]] = i;

    auto get = [&](int i)->std::string {
      return (i>=0 && i<(int)cells.size()) ? cells[i] : std::string();
    };

    int iK     = col(idx, {"k","strike"});
    int iT     = col(idx, {"t","maturity","expiry","tenor"});
    int iType  = col(idx, {"type","cp","callput"});

    int iBid   = col(idx, {"price_bid","bid","bid_price"});
    int iAsk   = col(idx, {"price_ask","ask","ask_price"});
    int iMid   = col(idx, {"mid","mid_price","price_mid"});

    int iIvBid = col(idx, {"iv_bid","sigma_bid","vol_bid"});
    int iIvAsk = col(idx, {"iv_ask","sigma_ask","vol_ask"});
    int iIvMid = col(idx, {"iv_mid","sigma_mid","vol_mid","implied_vol"});

    QuoteRow row;
    row.K       = parse_double(get(iK));
    row.T       = parse_double(get(iT));
    row.is_call = parse_is_call(get(iType));

    row.bid       = parse_double(get(iBid));
    row.ask       = parse_double(get(iAsk));
    row.mid_price = parse_double(get(iMid));

    row.iv_bid = parse_double(get(iIvBid));
    row.iv_ask = parse_double(get(iIvAsk));
    row.iv_mid = parse_double(get(iIvMid));

    // ---- Normalisations ----
    // mid prix si possible
    if (is_nan(row.mid_price) && is_finite(row.bid) && is_finite(row.ask) && row.ask >= row.bid) {
      row.mid_price = 0.5*(row.bid + row.ask);
    }
    // mid IV idem
    if (is_nan(row.iv_mid) && is_finite(row.iv_bid) && is_finite(row.iv_ask) && row.iv_ask >= row.iv_bid) {
      row.iv_mid = 0.5*(row.iv_bid + row.iv_ask);
    }

    // ---- Filtres ----
    bool bad = false;
    std::string why;

    if (!std::isfinite(row.K) || row.K <= 0.0) { bad = true; why = "K<=0 ou invalide"; }
    if (!bad && (!std::isfinite(row.T) || row.T <= 0.0)) { bad = true; why = "T<=0 ou invalide"; }

    // 1) prioriser les incohérences de spread (plus explicite)
    if (!bad && invalid_price_spread(row.bid, row.ask)) { bad = true; why = "spread prix aberrant"; }
    if (!bad && invalid_iv_spread(row.iv_bid, row.iv_ask)) { bad = true; why = "spread IV aberrant"; }

    // 2) ensuite seulement : “aucun mid exploitable”
    if (!bad && (std::isnan(row.mid_price) && std::isnan(row.iv_mid))) {
    bad = true; why = "ni prix mid ni IV mid";
    }

    if (bad) {
    if (num_ignored) (*num_ignored)++;
    if (warnings) warnings->push_back("Ligne ignorée: " + why);
    continue;
    }


    out.push_back(row);
  }

  return out;
}

} // namespace rw::io
