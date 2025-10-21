#pragma once
#include <vector>
#include <optional>
#include <cstddef>

namespace rw { namespace io { struct QuoteRow; } }

namespace rw::market {

// Surface d'IV basée sur w = sigma^2 * T, x = ln(K/S0)
class IvSurface {
public:
  IvSurface(double S0, double r, double q);

  // Ajout d'une quote (si iv_mid absent, on tentera price->IV lors du build)
  void add_quote(double K, double T, bool is_call,
                 std::optional<double> mid_price,
                 std::optional<double> iv_mid);

  // Ajout en masse depuis read_market_csv(...)
  void add_quotes_from_rows(const std::vector<rw::io::QuoteRow>& rows);

  // Construit les splines par maturité (PCHIP) + index temporel
  void build();

  // Evaluation IV(K,T) = sqrt( w(x,T) / T ), T>0. Retourne NaN si hors domaine.
  double iv(double K, double T) const;

  // Statistiques / debug
  std::size_t maturities() const;
  std::size_t quotes_kept() const { return kept_; }
  std::size_t quotes_skipped() const { return skipped_; }

private:
  struct Pending {
    double K;
    double T;
    bool   is_call;
    bool   has_price;
    double price;
    bool   has_iv;
    double iv;
  };

  struct Slice {
    double T{0.0};
    std::vector<double> x;   // ln(K/S0) croissant
    std::vector<double> w;   // variance implicite w = sigma^2 * T
    std::vector<double> d;   // dérivées PCHIP en x
    bool built{false};

    void build();            // calcule d (slopes) via PCHIP (Fritsch–Carlson)
    double eval(double xq) const; // evalue w(xq) avec cubic Hermite; extrapolation = clamp
  };

  double S0_;
  double r_;
  double q_;

  std::vector<Pending> pending_;   // quotes brutes ajoutées
  std::vector<Slice>   slices_;    // slices construites par maturité

  std::size_t kept_{0};
  std::size_t skipped_{0};

  static double clamp(double v, double lo, double hi);
};

} // namespace rw::market
