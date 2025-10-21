#pragma once
#include <cstddef>

namespace rw::market { struct MarketSurface; }
namespace rw { namespace smile  { struct SmileSurface; } }

namespace rw::calib {

struct CalibReport {
  double rmse_price{0.0};
  double rmse_iv{0.0};     // calculée si possible (quotes avec IV ou inversion)
  std::size_t n{0};        // nb de points utilisés
  bool converged{false};
  int iters{0};
};

// Calibre une seule sigma globale (BS homogène) sur toutes les quotes "priced" et arbitrage-OK.
// Retourne le rapport; écrit la sigma optimale dans sigma_opt.
// Bornes par défaut: [1e-6, 5.0], initialisation auto si sigma0 <= 0 (médiane des IV dispo sinon 0.2).
CalibReport fit_bs_global_sigma(const rw::market::MarketSurface& mkt,
                                double& sigma_opt,
                                double sigma0 = 0.0,
                                double sigma_min = 1e-6,
                                double sigma_max = 5.0,
                                int max_iters = 100,
                                double tol_step = 1e-10);
rw::smile::SmileSurface build_smile_surface(const rw::market::MarketSurface& mkt,
                                            double tolT = 1e-6);

} // namespace rw::calib
