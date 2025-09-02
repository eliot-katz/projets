#include <rw/core/stats.hpp>

#include <cmath>    // std::sqrt
#include <limits>   // std::numeric_limits

namespace rw {
namespace core {

// --- RunningStats -----------------------------------------------------------

RunningStats::RunningStats() noexcept = default;

void RunningStats::add(double x) noexcept {
  // Algorithme de Welford (stable numériquement, une passe)
  n_ += 1;
  const double delta  = x - mean_;
  mean_ += delta / static_cast<double>(n_);
  const double delta2 = x - mean_;
  m2_   += delta * delta2;
}

std::size_t RunningStats::count() const noexcept {
  return n_;
}

double RunningStats::mean() const noexcept {
  return mean_;
}

double RunningStats::variance() const noexcept {
  if (n_ < 2) {
    return std::numeric_limits<double>::quiet_NaN(); // politique: indéfini si n<2
  }
  return m2_ / static_cast<double>(n_ - 1); // variance d'échantillon
}

double RunningStats::std_error() const noexcept {
  if (n_ == 0) {
    return std::numeric_limits<double>::quiet_NaN(); // pas d'observations
  }
  const double var = variance(); // renvoie NaN si n<2, ce qui propage NaN ici si n==1
  return std::sqrt(var / static_cast<double>(n_));
}

// --- Confidence interval 95% ------------------------------------------------

[[nodiscard]] ConfidenceInterval confidence_interval_95(double mean,
                                                        double std_error,
                                                        std::size_t /*n*/) noexcept {
  // z pour 95% bilatéral sous normale
  static constexpr double Z95 = 1.959963984540054; 
  const double half = Z95 * std_error;
  return { mean - half, mean + half };
}

} // namespace core
} // namespace rw
