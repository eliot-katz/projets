#pragma once
/**
 * @file normal_rng.hpp
 * @brief Générateur de N(0,1) encapsulant uniquement la graine côté API.
 *
 * # Reproductibilité
 * Deux instances construites avec la même graine produisent la même séquence
 * de tirages. La copie/assignation copie la graine ; l'état interne est
 * reconstruit à partir de cette graine afin de préserver une ABI propre.
 *
 * # Concurrence
 * Utiliser **un RNG par thread**. Cette classe n'est pas thread-safe entre
 * threads pour une même instance.
 *
 * # Tests à prévoir
 * - Sur 1e6 tirages : moyenne ≈ 0 (|m| < ~3e-3) et variance ≈ 1 (|v-1| < ~3e-3).
 * - Reproductibilité : même graine ⇒ même séquence.
 * - Copie/assignation : la copie produit la même séquence que l’original
 *   reconstruit à partir de la graine.
 */

#include <cstdint>   // std::uint64_t
#include <cstddef>   // std::size_t (pour les tailles de buffers)
#include <vector>    // autorisé pour la sortie "bloc" si nécessaire (non exposé ici)
#include <memory>

namespace rw {
namespace core {

/**
 * @brief Générateur de lois normales standard N(0,1).
 *
 * Implémentation cachée (PIMPL) afin de ne pas exposer <random> dans l’API.
 * L’ABI reste stable : seule la graine est visible dans l’interface.
 */
class NormalRng {
public:
  /// @brief Construit avec une graine par défaut documentée.
  /// @details Par défaut, on utilise une graine fixe (implémentation) pour
  /// garantir la reproductibilité par défaut. Voir documentation d’implémentation.
  NormalRng();

  /// @brief Construit avec une graine explicite.
  explicit NormalRng(std::uint64_t seed);

  /// @brief Un tirage N(0,1).
  /// @return Un double ~ N(0,1).
  double sample() noexcept;

  /// @brief Remplit un buffer de n tirages N(0,1).
  /// @param out pointeur vers un buffer de taille >= n (non nul si n>0).
  /// @param n   nombre d’échantillons à produire.
  void sample_block(double* out, std::size_t n);

  /// @brief Graine utilisée pour (re)construire l’état interne.
  std::uint64_t seed() const noexcept;

  // Sémantique de copie/assignation : on copie la graine uniquement.
  NormalRng(const NormalRng&);
  NormalRng& operator=(const NormalRng&);

  // Déplacement par défaut.
  NormalRng(NormalRng&&) noexcept;
  NormalRng& operator=(NormalRng&&) noexcept;

  ~NormalRng() noexcept;

private:
  // PIMPL pour masquer <random> et garder l’ABI stable.
  struct Impl;
  std::unique_ptr<Impl> pimpl_;              // pointeur opaque vers l’implémentation
  std::uint64_t seed_;       // graine logique exposée via l’API
};

} // namespace core
} // namespace rw
