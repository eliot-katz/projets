#include <rw/core/normal_rng.hpp>

#include <random>   // std::mt19937_64, std::normal_distribution
#include <memory>   // std::make_unique
#include <cstddef>  // std::size_t (pour la boucle)

namespace rw {
namespace core {

namespace {
// Graine par défaut (documentée dans le header si tu veux) :
// constante liée au "golden ratio" de Knuth (bonne dispersion).
constexpr std::uint64_t DEFAULT_SEED = 0x9E3779B97F4A7C15ULL;
} // anonymous namespace

// --- PIMPL -------------------------------------------------------------------

struct NormalRng::Impl {
  std::mt19937_64 eng;
  std::normal_distribution<double> nd;

  Impl() : eng(), nd(0.0, 1.0) {}
};

// --- Ctors / Dtors -----------------------------------------------------------

NormalRng::NormalRng()
    : pimpl_(std::make_unique<Impl>()), seed_(DEFAULT_SEED) {
  pimpl_->eng.seed(seed_);
}

NormalRng::NormalRng(std::uint64_t seed)
    : pimpl_(std::make_unique<Impl>()), seed_(seed) {
  pimpl_->eng.seed(seed_);
}

// Copie "seed-only" : on reconstruit l'état à partir de la graine.
NormalRng::NormalRng(const NormalRng& other)
    : pimpl_(std::make_unique<Impl>()), seed_(other.seed_) {
  pimpl_->eng.seed(seed_);
}

NormalRng& NormalRng::operator=(const NormalRng& other) {
  if (this != &other) {
    seed_ = other.seed_;
    pimpl_ = std::make_unique<Impl>();
    pimpl_->eng.seed(seed_);
  }
  return *this;
}

// Déplacement : transfert du pimpl et de la graine.
NormalRng::NormalRng(NormalRng&& other) noexcept
    : pimpl_(std::move(other.pimpl_)), seed_(other.seed_) {}

NormalRng& NormalRng::operator=(NormalRng&& other) noexcept {
  if (this != &other) {
    pimpl_ = std::move(other.pimpl_);
    seed_  = other.seed_;
  }
  return *this;
}

NormalRng::~NormalRng() noexcept = default;

// --- API ---------------------------------------------------------------------

std::uint64_t NormalRng::seed() const noexcept {
  return seed_;
}

double NormalRng::sample() noexcept {
  return pimpl_->nd(pimpl_->eng);
}

void NormalRng::sample_block(double* out, std::size_t n) {
  if (n == 0) return;
  // Précondition (doc) : out != nullptr si n > 0.
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = pimpl_->nd(pimpl_->eng);
  }
}

} // namespace core
} // namespace rw
