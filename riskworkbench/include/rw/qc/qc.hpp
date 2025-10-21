#pragma once
#include <vector>

namespace rw { namespace market { struct MarketSurface; } }
namespace rw { namespace smile  { struct SmileSurface; } }

namespace rw::qc {

struct ErrorRow {
  double K{0.0}, T{0.0};
  bool   is_call{true};
  double price_mkt{0.0};
  double price_fit{0.0};
  double err_price{0.0};   // price_fit - price_mkt
  double iv_mkt{0.0};
  double iv_fit{0.0};
  double err_iv{0.0};      // iv_fit - iv_mkt
};

// Recalcule les prix via BS analytique à partir de iv_at(T,K) de la SmileSurface,
// construit un tableau de résidus (prix & IV). Ignore les lignes sans info exploitable.
std::vector<ErrorRow>
compare_fit(const rw::market::MarketSurface& mkt,
            const rw::smile::SmileSurface& smile);

} // namespace rw::qc
