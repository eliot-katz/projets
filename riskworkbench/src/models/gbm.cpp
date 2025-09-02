#include <rw/models/gbm.hpp>

#include <cmath>       // std::exp, std::sqrt
#include <stdexcept>   // std::invalid_argument

namespace rw {
namespace models {

Gbm::Gbm(GbmParams p) : params_(p) {
  if (p.sigma < 0.0) {
    throw std::invalid_argument("Gbm: sigma must be >= 0");
  }
}

double Gbm::sample_ST(double S0, double T, rw::core::NormalRng& rng) const {
  // Préconditions usuelles : S0 > 0, T >= 0. Cas limite :
  if (T == 0.0) {
    return S0;
  }
  const double sqrtT = std::sqrt(T);
  const double mu  = (params_.r - params_.q - 0.5 * params_.sigma * params_.sigma) * T;
  const double vol = params_.sigma * sqrtT;
  const double Z   = rng.sample();
  return S0 * std::exp(mu + vol * Z);
}

void Gbm::step_inplace(double& S, double dt, rw::core::NormalRng& rng) const {
  if (dt == 0.0) {
    return; // no-op
  }
  const double sqrt_dt = std::sqrt(dt);
  const double mu  = (params_.r - params_.q - 0.5 * params_.sigma * params_.sigma) * dt;
  const double vol = params_.sigma * sqrt_dt;
  const double Z   = rng.sample();
  S *= std::exp(mu + vol * Z);
}

} // namespace models
} // namespace rw
