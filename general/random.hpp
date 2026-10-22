#pragma once

#include "elff/general/error.hpp"
#include <array>
#include <cmath>
#include <cstdint>

namespace ELFF {
/**
 * @brief Stateless counter-based random stream for reproducible model setup.
 */
class CounterRandom {
public:
  static constexpr std::uint64_t default_seed = 1;

  explicit constexpr CounterRandom(std::uint64_t seed = default_seed,
                                  std::uint64_t stream = 0)
      : seed_(seed), stream_(stream) {}
  constexpr std::uint64_t seed() const { return seed_; }
  constexpr std::uint64_t stream_id() const { return stream_; }

  /** @brief Original SplitMix64-based integer mixing operation. */
  static constexpr std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }
  /** @brief Combine stable identifiers, e.g. model ID and timestep index. */
  static constexpr std::uint64_t fold(std::uint64_t a, std::uint64_t b) {
    return mix64(a ^ mix64(b));
  }
  /** @brief Derive a child stream without consuming draws or changing this stream.
   * The child retains the master seed and uses fold(stream_id(), id).
   */
  constexpr CounterRandom derive(std::uint64_t id) const {
    return CounterRandom(seed_, fold(stream_, id));
  }
  /** @brief Original deterministic seed/stream/index mapping. */
  constexpr std::uint64_t u64(std::uint64_t index) const {
    return fold(fold(seed_, stream_), index);
  }
  /** @brief Uniform deviate in \f([0,1)\f), using the upper 53 random bits. */
  double uniform01(std::uint64_t index) const {
    return static_cast<double>(u64(index) >> 11) * 0x1.0p-53;
  }
  /** @brief Uniform deviate in \f((0,1)\f); replace zero by \f(2^{-53}\f). */
  double uniform01_open(std::uint64_t index) const {
    const double u = uniform01(index);
    return u > 0. ? u : 0x1.0p-53;
  }
  /** @brief Standard normal using the original Box-Muller counters and transform.
   * \f[z=\sqrt{-2\log u_1}\cos(2\pi u_2),\quad
   * u_1=U_{\mathrm{open}}(2i),\ u_2=U(2i+1).\f]
   */
  double normal(std::uint64_t index) const {
    const double u1 = uniform01_open(2 * index);
    const double u2 = uniform01(2 * index + 1);
    constexpr double pi = 3.14159265358979323846;
    return std::sqrt(-2. * std::log(u1)) * std::cos(2. * pi * u2);
  }
  /** @brief Normal vector with counters \f(3i+d\f); inactive coordinates are zero. */
  std::array<double, 3> gaussian_vector(std::uint64_t index, unsigned dimensions = 3) const {
    ELFF_ASSERT(dimensions >= 1 && dimensions <= 3, "Random vector dimensions must be 1, 2 or 3.");
    std::array<double, 3> result{};
    for (unsigned d = 0; d < dimensions; ++d) result[d] = normal(3 * index + d);
    return result;
  }
  /** @brief Uniform box position with the original coordinate counter mapping.
   * Active coordinates are origin[d] + extent[d] * U(3i+d); inactive ones are
   * zero, matching ParticleRandom. Supply physical domain extents explicitly.
   */
  std::array<double, 3> domain_position(std::uint64_t index,
      const std::array<double, 3>& origin, const std::array<double, 3>& extent,
      unsigned dimensions = 3) const {
    ELFF_ASSERT(dimensions >= 1 && dimensions <= 3, "Random position dimensions must be 1, 2 or 3.");
    std::array<double, 3> result{};
    for (unsigned d = 0; d < dimensions; ++d) {
      ELFF_ASSERT(std::isfinite(origin[d]) && std::isfinite(extent[d]) && extent[d] >= 0 &&
                  std::isfinite(origin[d] + extent[d]), "Invalid random position domain.");
      result[d] = origin[d] + extent[d] * uniform01(3 * index + d);
    }
    return result;
  }
private:
  std::uint64_t seed_; /**< @brief Master seed, unchanged by stream derivation. */
  std::uint64_t stream_; /**< @brief Stable stream identifier. */
};
} // namespace ELFF
