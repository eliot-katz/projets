#include "rw/io/market_csv.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm> // any_of
using namespace std;

int main(int argc, char** argv) {
  const string path = (argc>1 ? argv[1] : "data/market_samples/sample_market.csv");

  size_t ignored = 0;
  vector<string> warnings;
  auto rows = rw::io::read_market_csv(path, &ignored, &warnings);

  cout << "Valid rows: " << rows.size() << "\n";
  cout << "Ignored rows: " << ignored << "\n";

  if (!rows.empty()) {
    const auto& r0 = rows.front();
    cout << "First valid: K=" << r0.K
         << " T=" << r0.T
         << " is_call=" << (r0.is_call?"C":"P")
         << " mid_price=" << r0.mid_price
         << " iv_mid=" << r0.iv_mid << "\n";
  }

  // Asserts “attendus”
  constexpr double EPS = 1e-12;
  assert(rows.size() == 6);
  assert(ignored == 4);
  assert(std::abs(rows.front().mid_price - 100.0) < EPS);
  assert(rows.front().is_call == true);

  // Vérifie qu’on a bien des messages pour K, T, spread prix, et “ni mid”
  auto has_warn = [&](const string& needle){
    return any_of(warnings.begin(), warnings.end(),
                  [&](const string& w){ return w.find(needle) != string::npos; });
  };
  assert(has_warn("K<=0 ou invalide"));
  assert(has_warn("T<=0 ou invalide"));
  assert(has_warn("spread prix aberrant"));
  assert(has_warn("ni prix mid ni IV mid"));

  for (auto& w: warnings) cerr << "[warn] " << w << "\n";
  return 0;
}
