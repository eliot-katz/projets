#include "rw/io/market_csv.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

int main(int argc, char** argv) {
  std::string path;
  bool show_warnings = false;
  bool require_meta  = false;

  for (int i=1;i<argc;++i) {
    std::string a = argv[i];
    if ((a=="-f" || a=="--file") && i+1<argc) { path = argv[++i]; }
    else if (a=="-w" || a=="--show-warnings") { show_warnings = true; }
    else if (a=="-m" || a=="--require-meta")  { require_meta  = true; }
    else if (a=="-h" || a=="--help") {
      std::cout << "Usage: market_csv_info -f <file.csv> [-w] [-m]\n";
      return 0;
    } else if (path.empty()) { path = a; }
  }
  if (path.empty()) {
    std::cerr << "Please provide a CSV path (-f <file.csv>).\n";
    return 2;
  }

  std::size_t ignored = 0;
  std::vector<std::string> warnings;
  auto rows = rw::io::read_market_csv(path, &ignored, &warnings);

  std::cout << "File: " << path << "\n";
  std::cout << "Valid rows: " << rows.size() << "\n";
  std::cout << "Ignored rows: " << ignored << "\n";

  if (require_meta) {
    auto meta = rw::io::read_market_meta_csv(path);
    if (!meta || !std::isfinite(meta->S0) || !std::isfinite(meta->r) || !std::isfinite(meta->q)) {
      std::cerr << "Missing market meta (underlying,S0,r,q).\n";
      return 3;
    }
    std::cout << "Meta: underlying=" << meta->underlying
              << " S0=" << meta->S0 << " r=" << meta->r << " q=" << meta->q << "\n";
  }

  if (show_warnings) {
    for (auto& w : warnings) std::cerr << "[warn] " << w << "\n";
  }
  return 0;
}
