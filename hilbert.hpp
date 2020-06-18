#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <cstdint>

namespace mve {

uint64_t hilbertToMorton3(const uint64_t hilbert, const unsigned iteration) noexcept;
uint64_t mortonToHilbert3(const uint64_t morton, const unsigned iteration) noexcept;

uint32_t fastMortonToHilbert3_32(const uint32_t morton, const uint32_t bits) noexcept;
uint32_t fastHilbertToMorton3_32(const uint32_t morton, const uint32_t bits) noexcept;

}  // namespace mve

#endif  // HILBERT_HPP
