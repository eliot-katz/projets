#pragma once
#include <cstddef>

namespace rw::bs {

struct IvResult {
  double sigma;     // solution (ou borne si non convergé)
  int    iters;     // itérations effectuées (Newton + éventuelle bisection)
  bool   converged; // true si tolérance atteinte
};

// Prix -> IV (CALL par défaut). Bornes: [1e-6, 5.0]
IvResult implied_vol(double S0, double K, double r, double q, double T, double price);

// Variante explicite CALL/PUT (PUT converti en CALL via parité)
IvResult implied_vol_cp(double S0, double K, double r, double q, double T, double price, bool is_call);

} // namespace rw::bs
