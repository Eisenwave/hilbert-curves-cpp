#include "hilbert.hpp"

#include <array>
#include <set>

namespace mve {

using permutation8 = std::array<uint8_t, 8>;

template <typename T>
using vec3 = std::array<T, 3>;
using vec3i = vec3<int>;
using vec3u = vec3<unsigned>;

// ======== USER CONFIGURABLE AREA =====================================================================================

static constexpr permutation8 HILBERT_PATTERN_A = {0, 1, 3, 2, 6, 7, 5, 4};
static constexpr permutation8 HILBERT_PATTERN_B = {0, 2, 6, 7, 3, 1, 5, 4};
static constexpr permutation8 HILBERT_PATTERN_C = {0, 1, 5, 4, 6, 2, 3, 7};

// currently only works with gray codes, pattern C is theoretically viable too but not a gray code
static constexpr size_t PATTERN = 0;

/**
 * @brief The hilbert curve base pattern. This must be a gray code.
 * This permutation maps each index of the hilbert curve onto a morton (bit-interleaved) position in 3d space.
 */
constexpr permutation8 RAW_HILBERT_TO_MORTON =
    PATTERN == 0 ? HILBERT_PATTERN_A : PATTERN == 1 ? HILBERT_PATTERN_B : HILBERT_PATTERN_C;

/**
 * @brief ROT_HILBERT_TO_MORTON the rotation of the hilbert permutation to the left or right.
 * Hilbert curves require a gray-code as a base pattern, but the order of the elements can be chosen at will.
 *
 * Change this value to rotate to the right (if positive) or to the left (if negative).
 * This constant exists purely for convenience, so that the user doesn't have to manually change RAW_HILBERT_TO_MORTON.
 */
constexpr int ROT_HILBERT_TO_MORTON = 0;

// Uncomment this if you don't want to recompute the permutations
// You will have to manually define them yourself in the #else clause.
//
// #define DONT_CONSTRUCT_PERMUTATIONS

// #define DO_PRINT_CONSTRUCTED_PERMUTATIONS
// =====================================================================================================================

template <typename T>
constexpr vec3<T> operator+(vec3<T> l, vec3<T> r)
{
    return {l[0] + r[0], l[1] + r[1], l[2] + r[2]};
}

template <typename T>
constexpr vec3<T> operator-(vec3<T> l, vec3<T> r)
{
    return {l[0] - r[0], l[1] - r[1], l[2] - r[2]};
}

template <typename T>
constexpr vec3<T> operator*(vec3<T> v, T s)
{
    return {v[0] * s, v[1] * s, v[2] * s};
}

template <typename T>
constexpr vec3<T> operator/(vec3<T> v, T s)
{
    return {v[0] / s, v[1] / s, v[2] / s};
}

template <typename T>
static constexpr T dot(vec3<T> x, vec3<T> y)
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

#ifndef DONT_CONSTRUCT_PERMUTATIONS
static constexpr vec3i matrixMultiply(vec3i rows[3], vec3i v)
{
    return {dot(rows[0], v), dot(rows[1], v), dot(rows[2], v)};
}

/**
 * @brief Compute the inverse p^-1 to a permutation p so that p^-1 * p = id.
 * @param perm the permutation to invert
 * @return the inverse permutation
 */
static constexpr permutation8 invertPerm(const permutation8 &perm)
{
    permutation8 inverse{};
    for (uint8_t val = 0; val < 8; ++val) {
        for (uint8_t i = 0; i < 8; ++i) {
            if (perm[i] == val) {
                inverse[val] = i;
                break;
            }
        }
    }
    return inverse;
}

/**
 * Unsigned modulo operation.
 * This will always return a positive result.
 * For instance, umod(-1, 16) = 15.
 * @param n the number
 * @param mod the modulus
 * @return n % m
 */
template <typename Int,
          typename Uint = std::make_unsigned_t<Int>,
          std::enable_if_t<(std::is_integral_v<Int> && std::is_unsigned_v<Uint>), int> = 0>
constexpr Uint umod(Int n, Uint mod) noexcept
{
    if constexpr (std::is_unsigned_v<Int>)
        return n % mod;
    else
        return static_cast<Uint>(n % static_cast<Int>(mod) + static_cast<Int>(mod)) % mod;
}

/**
 * @brief Rotates a permutation left or right by a given offset.
 * If the offset is positive, the permutation will be rotated to the right.
 * @param perm the permutation to rotate
 * @return the rotated permutation
 */
static constexpr permutation8 rotatePerm(const permutation8 &perm, int offset)
{
    permutation8 result{};
    for (unsigned source = 0; source < 8; ++source) {
        unsigned target = umod(static_cast<int>(source) + offset, 8u);
        result[target] = perm[source];
    }
    return result;
}

constexpr permutation8 hilbertToMortonPerm = rotatePerm(RAW_HILBERT_TO_MORTON, ROT_HILBERT_TO_MORTON);
constexpr permutation8 mortonToHilbertPerm = invertPerm(hilbertToMortonPerm);

static_assert(mortonToHilbertPerm[hilbertToMortonPerm[0]] == 0);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[1]] == 1);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[2]] == 2);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[3]] == 3);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[4]] == 4);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[5]] == 5);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[6]] == 6);
static_assert(mortonToHilbertPerm[hilbertToMortonPerm[7]] == 7);

static constexpr permutation8 matrixToPerm(vec3i rows[3])
{
    permutation8 result{};
    for (uint8_t i = 0; i < 8; ++i) {
        const vec3i vec = vec3i{(i >> 2) & 1, (i >> 1) & 1, (i >> 0) & 1} * 2 - vec3i{1, 1, 1};

        vec3i rotated = (matrixMultiply(rows, vec) + vec3i{1, 1, 1}) / 2;

        result[i] = static_cast<uint8_t>((rotated[0] << 2) | (rotated[1] << 1) | (rotated[2] << 0));
    }

    return result;
}

#ifdef MVE_DEBUG
static constexpr bool isPermutation(const permutation8 &perm)
{
    for (size_t i = 0; i < 8; ++i) {
        if (perm[i] > 7) {
            return false;
        }
        for (size_t j = 0; j < 8; ++j) {
            if (i != j && perm[i] == perm[j]) {
                return false;
            }
        }
    }
    return true;
}
#endif

static constexpr permutation8 makeRotationPermutation(size_t rx, size_t ry, size_t rz)
{
    constexpr int sin[4]{0, 1, 0, -1};
    constexpr int cos[4]{1, 0, -1, 0};

    int sx = sin[rx];
    int sy = sin[ry];
    int sz = sin[rz];

    int cx = cos[rx];
    int cy = cos[ry];
    int cz = cos[rz];

    // https://wikimedia.org/api/rest_v1/media/math/render/svg/ecef3eab0b5f54b66c7107f3ceeeb319a2e0e756
    vec3i rotRows[3]{{cx * cy, sx * sz - cx * cz * sy, cz * sx + cx * sy * sz},
                     {sy, cy * cz, -cy * sz},
                     rotRows[2] = {-cy * sx, cx * sz + cz * sx * sy, cx * cz - sx * sy * sz}};

    return matrixToPerm(rotRows);
}

static std::array<permutation8, 24> makeAllRotationPermutations()
{
    std::set<permutation8> setOfRotations;

    for (size_t rx = 0; rx < 4; ++rx) {
        for (size_t ry = 0; ry < 4; ++ry) {
            for (size_t rz = 0; rz < 4; ++rz) {
                auto perm = makeRotationPermutation(rx, ry, rz);
                setOfRotations.emplace(perm);
            }
        }
    }

    std::array<permutation8, 24> result;
    size_t index = 0;
    for (auto pos = setOfRotations.begin(); pos != setOfRotations.end(); ++pos) {
        result[index++] = std::move(*pos);
    }

    return result;
}

static std::array<permutation8, 24> allRotationPerms = makeAllRotationPermutations();

static constexpr bool areDirectNeighbors(vec3u a, vec3u b)
{
    vec3i diff = vec3i{static_cast<int>(a[0]) - static_cast<int>(b[0]),
                       static_cast<int>(a[1]) - static_cast<int>(b[1]),
                       static_cast<int>(a[2]) - static_cast<int>(b[2])};
    return dot(diff, diff) == 1;
}

static constexpr std::array<vec3u, 8> makePositions(const permutation8 &code)
{
    std::array<vec3u, 8> result{};
    for (uint8_t i = 0; i < 8; ++i) {
        uint8_t val = code[i];

        unsigned x = (val >> 2) & 1;
        unsigned y = (val >> 1) & 1;
        unsigned z = (val >> 0) & 1;
        result[i] = {x, y, z};
    }
    return result;
}

constexpr permutation8 identityPerm = {0, 1, 2, 3, 4, 5, 6, 7};

constexpr std::array<vec3u, 8> mortonPositions = makePositions(identityPerm);
constexpr std::array<vec3u, 8> hilbertToMortonPositions = makePositions(hilbertToMortonPerm);

template <bool INVERSE>
static constexpr bool makeHilbertToMortonPermutations_impl(permutation8 out[8], size_t i, vec3u prev_end_global)
{
    constexpr size_t last_index = 7;
    constexpr vec3u begin_global = hilbertToMortonPositions[0] * 2u + hilbertToMortonPositions[0];
    constexpr vec3u end_global = hilbertToMortonPositions[last_index] * 2u + hilbertToMortonPositions[last_index];

    const vec3u curr_global = hilbertToMortonPositions[i] * 2u;

    const auto emitResult = [&out, i](size_t permIndex) {
        const size_t outIndex = hilbertToMortonPerm[i];
        if constexpr (INVERSE) {
            out[outIndex] = invertPerm(allRotationPerms[permIndex]);
        }
        else {
            out[outIndex] = allRotationPerms[permIndex];
        }
    };

    // we try out all 24 possible rotation-permutations until we find one which fits
    // in this context "fits" means that it seamlessly connects to the previous piece
    // the first piece must align with begin_global
    // the last piece must align with end_global
    for (size_t perm_index = 0; perm_index < 24; ++perm_index) {
        uint8_t hilbert_begin = allRotationPerms[perm_index][hilbertToMortonPerm[0]];
        vec3u curr_begin_local = mortonPositions[hilbert_begin];
        vec3u curr_begin_global = curr_begin_local + curr_global;

        if (i == 0) {
            // we know the global begin and this is what we test at the first iteration
            if (curr_begin_global != begin_global) {
                continue;
            }
        }
        else {
            // the begin of our current 2x2x2 block has to seamlessly connect with the previous block
            if (not areDirectNeighbors(curr_begin_global, prev_end_global)) {
                continue;
            }
        }

        uint8_t hilbert_end = allRotationPerms[perm_index][hilbertToMortonPerm[last_index]];
        vec3u curr_end_local = mortonPositions[hilbert_end];
        vec3u curr_end_global = curr_end_local + curr_global;

        // unless we reached the end, we must try to connect the next piece
        if (i != last_index) {
            if (makeHilbertToMortonPermutations_impl<INVERSE>(out, i + 1, curr_end_global)) {
                // yay, we succeeded, now we fill up the result and forward true
                emitResult(perm_index);
                return true;
            }
            else {
                continue;
            }
        }

        // this test will be performed for the final piece
        if (curr_end_global != end_global) {
            continue;
        }

        // yay, we found a complete permutation array from begin to end, now fill up the result
        emitResult(perm_index);
        return true;
    }

    return false;
}

static std::array<permutation8, 8> makeHilbertPermutations()
{
    std::array<permutation8, 8> result{};
    makeHilbertToMortonPermutations_impl<false>(result.data(), 0, vec3u{0, 0, 0});

    return result;
}

static constexpr std::array<permutation8, 8> invertPerms(const std::array<permutation8, 8> &perms)
{
    std::array<permutation8, 8> result{};
    for (size_t i = 0; i < 8; ++i) {
        result[i] = invertPerm(perms[i]);
    }
    return result;
}

static const std::array<permutation8, 8> hilbertToMortonPermutations = makeHilbertPermutations();
static const std::array<permutation8, 8> mortonToHilbertPermutations = invertPerms(hilbertToMortonPermutations);

#ifdef DO_PRINT_CONSTRUCTED_PERMUTATIONS
static void printPermutations(const char *name, const permutation8 perms[8])
{
    std::cerr << "static constexpr permutation8 ";
    std::cerr << name;
    std::cerr << "[8] = {\n";
    for (size_t i = 0; i < 8; ++i) {
        std::string row = "    {";
        for (size_t j = 0; j < 8; ++j) {
            row += std::to_string(perms[i][j]);
            if (j != 7) {
                row += ", ";
            }
        }
        row += "}\n";
        std::cerr << row;
    }
    std::cerr << "}\n";
}

static int doPrintPermutations()
{
    printPermutations("hilbertToMortonPermutations", hilbertToMortonPermutations.data());
    printPermutations("mortonToHilbertPermutations", mortonToHilbertPermutations.data());
    return 0;
}

static const int staticMain = doPrintPermutations();
#endif

#else
static constexpr permutation8 hilbertToMortonPermutations[8] = {
    {0, 4, 1, 5, 2, 6, 3, 7},
    {0, 2, 4, 6, 1, 3, 5, 7},
    {3, 7, 2, 6, 1, 5, 0, 4},
    {6, 4, 2, 0, 7, 5, 3, 1},
    {3, 1, 7, 5, 2, 0, 6, 4},
    {0, 1, 2, 3, 4, 5, 6, 7},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {6, 4, 2, 0, 7, 5, 3, 1}
};
static constexpr permutation8 mortonToHilbertPermutations[8] = {
    {0, 2, 4, 6, 1, 3, 5, 7},
    {0, 4, 1, 5, 2, 6, 3, 7},
    {6, 4, 2, 0, 7, 5, 3, 1},
    {3, 7, 2, 6, 1, 5, 0, 4},
    {5, 1, 4, 0, 7, 3, 6, 2},
    {0, 1, 2, 3, 4, 5, 6, 7},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {3, 7, 2, 6, 1, 5, 0, 4},
};
#endif

uint64_t hilbertToMorton3(const uint64_t hilbert, const unsigned iteration) noexcept
{
    const uint8_t leastHilbertDigit = hilbert & 0b111;
    const uint8_t leastMortonDigit = hilbertToMortonPerm[leastHilbertDigit];
    uint64_t result = (hilbert ^ uint64_t{leastHilbertDigit}) | uint64_t{leastMortonDigit};

    for (unsigned i = 1; i <= iteration; ++i) {
        const unsigned shift = i * 3;
        const uint8_t hilbertDigit = (result >> shift) & 0b111;
        const uint8_t mortonDigit = hilbertToMortonPerm[hilbertDigit];
        result ^= uint64_t{hilbertDigit} << shift;  // clear the original digit
        result |= uint64_t{mortonDigit} << shift;   // fill in the new digit

        const permutation8 perm = hilbertToMortonPermutations[mortonDigit];

        for (unsigned j = i - 1; j < i; --j) {
            const unsigned childShift = j * 3;
            const uint8_t mortonChildDigit = (result >> childShift) & 0b111;
            const uint8_t permutedChildDigit = perm[mortonChildDigit];
            result ^= uint64_t{mortonChildDigit} << childShift;    // clear the original digit
            result |= uint64_t{permutedChildDigit} << childShift;  // fill in the new digit
        }
    }

    return result;
}

uint64_t mortonToHilbert3(const uint64_t morton, const unsigned iteration) noexcept
{
    uint64_t result = morton;

    for (unsigned i = iteration; i != 0; --i) {
        const unsigned shift = i * 3;
        const uint8_t mortonDigit = (result >> shift) & 0b111;
        const uint8_t hilbertDigit = mortonToHilbertPerm[mortonDigit];
        result ^= uint64_t{mortonDigit} << shift;   // clear the original digit
        result |= uint64_t{hilbertDigit} << shift;  // fill in the new digit

        const permutation8 perm = mortonToHilbertPermutations[mortonDigit];

        for (unsigned j = i - 1; j < i; --j) {
            const unsigned childShift = j * 3;
            const uint8_t mortonChildDigit = (result >> childShift) & 0b111;
            const uint8_t permutedChildDigit = perm[mortonChildDigit];
            result ^= uint64_t{mortonChildDigit} << childShift;    // clear the original digit
            result |= uint64_t{permutedChildDigit} << childShift;  // fill in the new digit
        }
    }

    const uint8_t leastMortonDigit = result & 0b111;
    const uint8_t leastHilbertDigit = mortonToHilbertPerm[leastMortonDigit];
    return (result ^ leastMortonDigit) | leastHilbertDigit;
}

uint32_t fastMortonToHilbert3_32(const uint32_t morton, const uint32_t bits) noexcept
{
    uint32_t hilbert = morton;
    if (bits > 1) {
        uint32_t block = ((bits * 3) - 3);
        uint32_t hcode = ((hilbert >> block) & 7);
        uint32_t mcode, shift, signs;
        shift = signs = 0;
        while (block) {
            block -= 3;
            hcode <<= 2;
            mcode = ((0x20212021 >> hcode) & 3);
            shift = ((0x48 >> (7 - shift - mcode)) & 3);
            signs = ((signs | (signs << 3)) >> mcode);
            signs = ((signs ^ (0x53560300 >> hcode)) & 7);
            mcode = ((hilbert >> block) & 7);
            hcode = mcode;
            hcode = (((hcode | (hcode << 3)) >> shift) & 7);
            hcode ^= signs;
            hilbert ^= ((mcode ^ hcode) << block);
        }
    }
    hilbert ^= ((hilbert >> 1) & 0x92492492);
    hilbert ^= ((hilbert & 0x92492492) >> 1);
    return (hilbert);
}

uint32_t fastHilbertToMorton3_32(const uint32_t hilbert, const uint32_t bits) noexcept
{
    uint32_t morton = hilbert;
    morton ^= ((morton & 0x92492492) >> 1);
    morton ^= ((morton >> 1) & 0x92492492);
    if (bits > 1) {
        uint32_t block = ((bits * 3) - 3);
        uint32_t hcode = ((morton >> block) & 7);
        uint32_t mcode, shift, signs;
        shift = signs = 0;
        while (block) {
            block -= 3;
            hcode <<= 2;
            mcode = ((0x20212021 >> hcode) & 3);
            shift = ((0x48 >> (4 - shift + mcode)) & 3);
            signs = ((signs | (signs << 3)) >> mcode);
            signs = ((signs ^ (0x53560300 >> hcode)) & 7);
            hcode = ((morton >> block) & 7);
            mcode = hcode;
            mcode ^= signs;
            mcode = (((mcode | (mcode << 3)) >> shift) & 7);
            morton ^= ((hcode ^ mcode) << block);
        }
    }
    return (morton);
}

}  // namespace mve
