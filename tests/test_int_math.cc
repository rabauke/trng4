// Copyright (c) 2000-2026, Heiko Bauke
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.

#include <array>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <cmath>
#if defined _MSC_VER && __cplusplus <= 201703
#include <ciso646>
#endif

#include <catch2/catch_all.hpp>

#include <trng/int_math.hpp>

TEST_CASE("integer math") {
  SECTION("pow2") {
    REQUIRE(trng::int_math::pow2(0) == 1);
    REQUIRE(trng::int_math::pow2(8) == 256);
    REQUIRE(trng::int_math::pow2(63ull) == 0x8000000000000000ull);
  }

  SECTION("log2_floor") {
    REQUIRE(trng::int_math::log2_floor(1) == 0);
    REQUIRE(trng::int_math::log2_floor(0x10000) == 16);
    REQUIRE(trng::int_math::log2_floor(0x1ffff) == 16);
  }

  SECTION("log2_ceil") {
    REQUIRE(trng::int_math::log2_ceil(1) == 0);
    REQUIRE(trng::int_math::log2_ceil(0x10000) == 16);
    REQUIRE(trng::int_math::log2_ceil(0x10001) == 17);
  }

  SECTION("ceil2") {
    REQUIRE(trng::int_math::ceil2(1) == 1);
    REQUIRE(trng::int_math::ceil2(0x10000) == 0x10000);
    REQUIRE(trng::int_math::ceil2(0x10001) == 0x20000);
  }

  SECTION("mask") {
    REQUIRE(trng::int_math::mask(0x10000) == 0x1ffff);
    REQUIRE(trng::int_math::mask(0x10001) == 0x1ffff);
    REQUIRE(trng::int_math::mask(trng::int32_t(0x7fffffff)) == 0x7fffffff);
    REQUIRE(trng::int_math::mask(trng::uint32_t(0xffffffffu)) == 0xffffffffu);
  }
}


TEST_CASE("linear algebra") {
  SECTION("matrix vector multiplication") {
    const int n{3};
    // clang-format off
    const trng::int32_t A[n * n]{
      1, 2, 3,
      2, 4, 6,
      3, 2, 5 };
    // clang-format on
    const trng::int32_t b[n]{2, 3, 4};
    const trng::int32_t m{7};
    trng::int32_t c[n];
    const trng::int32_t c_exact[n]{6, 5, 4};
    trng::int_math::matrix_vec_mult<n>(A, b, c, m);
    for (int i{0}; i < n; ++i)
      REQUIRE(c[i] == c_exact[i]);
  }

  SECTION("matrix matrix multiplication") {
    const int n{3};
    // clang-format off
    const trng::int32_t A[n * n]{
      1, 2, 3,
      2, 4, 6,
      3, 2, 5 };
    const trng::int32_t B[n * n]{
      2, 2, 3,
      2, 3, 6,
      3, 2, 4 };
    const trng::int32_t C_exact[n * n] {
      1, 0, 6,
      2, 0, 5,
      4, 1, 6 };
    // clang-format on
    const trng::int32_t m{7};
    trng::int32_t C[n * n];
    trng::int_math::matrix_mult<3>(A, B, C, m);
    for (int i{0}; i < n * n; ++i)
      REQUIRE(C[i] == C_exact[i]);
  }

  SECTION("Gaussian elimination") {
    const int n{3};
    // clang-format off
    trng::int32_t A[n * n]{
      1, 2, 3,
      2, 1, 6,
      3, 2, 5 };
    // clang-format on
    trng::int32_t b[n]{2, 4, 3};
    const trng::int32_t m{7};
    const trng::int32_t x_exact[n]{5, 0, 6};
    trng::int_math::gauss<n>(A, b, m);
    for (int i{0}; i < n; ++i)
      REQUIRE(b[i] == x_exact[i]);
  }

  SECTION("Gaussian elimination with singular matrix") {
    const int n{3};
    // clang-format off
    trng::int32_t A[n * n]{
      1, 2, 3,
      2, 4, 6,
      3, 2, 5 };
    const trng::int32_t A_orig[n * n]{
      1, 2, 3,
      2, 4, 6,
      3, 2, 5 };
    // clang-format on
    trng::int32_t b[n]{2, 4, 3};
    const trng::int32_t b_orig[n]{2, 4, 3};
    const trng::int32_t m{7};
    // solution is not unique because the matrix is singular, gaussian elimination picks one
    // solution if one exists
    trng::int_math::gauss<n>(A, b, m);
    trng::int32_t c[n]{0, 0, 0};
    trng::int_math::matrix_vec_mult<3>(A_orig, b, c, m);
    for (int i{0}; i < n; ++i)
      REQUIRE(c[i] == b_orig[i]);
  }

  SECTION("Gaussian elimination with singular matrix and no solution") {
    const int n{3};
    // clang-format off
    trng::int32_t A[n * n]{
      1, 2, 3,
      2, 4, 6,
      3, 2, 5 };
    // clang-format on
    trng::int32_t b[n]{2, 1, 3};
    const trng::int32_t m{7};
    // solution does not exist because the matrix is singular and the choice of the right-hand
    // side
    REQUIRE_THROWS_AS(trng::int_math::gauss<n>(A, b, m), std::exception);
  }
}


TEST_CASE("modulo inverse") {
  SECTION("prime modulus") {
    const long m{104729};  // must be prime
    for (long a{1}; a < m; ++a) {
      int b{trng::int_math::modulo_inverse(a, m)};
      REQUIRE((static_cast<long long>(a) * static_cast<long long>(b)) % m == 1);
    }
    REQUIRE_THROWS_AS(trng::int_math::modulo_inverse(0, m), std::exception);
  }

  SECTION("power-of-2 modulus") {
    const long m{1024 * 1024};  // must be power of 2
    for (long a{1}; a < m; a += 2) {
      int b{trng::int_math::modulo_inverse(a, m)};
      REQUIRE((static_cast<long long>(a) * static_cast<long long>(b)) % m == 1);
    }
    REQUIRE_THROWS_AS(trng::int_math::modulo_inverse(100, m), std::exception);
  }
}
