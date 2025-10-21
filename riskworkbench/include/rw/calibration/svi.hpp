#pragma once
#include <vector>
#include <cstddef>

namespace rw { namespace market { struct MarketSurface; } }

namespace rw::svi {

// Paramétrisation raw-SVI: w(k) = a + b * ( rho*(k - m) + sqrt((k - m)^2 + sigma^2 ) )
struct Params {
  double a{0.01};     // niveau
  double b{0.10};     // pente ailes (>=0)
  double rho{0.0};    // skew (|rho|<1)
  double m{0.0};      // translation (k-shift)
  double sigma{0.20}; // largeur noyau (>0)
};

struct FitReport {
  bool converged{false};
  int  iters{0};
  std::size_t n{0};
  double rmse_w{0.0};      // RMSE sur w
  double rmse_iv{0.0};     // RMSE sur sigma (si on a des IV de référence)
};

// Évalue w_SVI(k; p). Renvoie NaN si paramètres invalides.
double w_svi(double k, const Params& p);

// Fit d’une “slice” (une maturité): on donne k = ln(K/S0), w = sigma^2 * T.
// Retourne le rapport et écrit les paramètres dans out.
// Réglages: itérations LM simples + bornes/clamp (b>=1e-8, sigma>=1e-6, |rho|<=0.999).
FitReport fit_slice(const std::vector<double>& k,
                    const std::vector<double>& w,
                    Params& out,
                    int max_iters = 100,
                    double tol_step = 1e-10);

// ---- Surface SVI (une set de params par maturité) ----
struct Slice {
  double T{0.0};
  Params p;
};

struct Surface {
  double S0{1.0}, r{0.0}, q{0.0};
  std::vector<Slice> slices; // trié par T croissant

  // iv(K,T) = sqrt( w_svi( ln(K/S0); p_T ) / T )
  double iv(double K, double T) const;
};

} // namespace rw::svi
