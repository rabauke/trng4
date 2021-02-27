// Copyright (c) 2000-2021, Heiko Bauke
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
#include <ciso646>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/int_math.hpp>


BOOST_AUTO_TEST_SUITE(test_suite_int_math)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_integer_math)

BOOST_AUTO_TEST_CASE(test_pow2) {
  BOOST_TEST(trng::int_math::pow2(0) == 1, "2^0 = 1");
  BOOST_TEST(trng::int_math::pow2(8) == 256, "2^8 = 256");
  BOOST_TEST(trng::int_math::pow2(63ull) == 0x8000000000000000ull,
             "2^63 = 9223372036854775808");
}

BOOST_AUTO_TEST_CASE(test_log2_floor) {
  BOOST_TEST(trng::int_math::log2_floor(1) == 0, "ld_floor 1 = 0");
  BOOST_TEST(trng::int_math::log2_floor(0x10000) == 16, "ld_floor 65536 = 16");
  BOOST_TEST(trng::int_math::log2_floor(0x1ffff) == 16, "ld_floor 131071 = 16");
}

BOOST_AUTO_TEST_CASE(test_log2_ceil) {
  BOOST_TEST(trng::int_math::log2_ceil(1) == 0, "ld_ceil 1 = 0");
  BOOST_TEST(trng::int_math::log2_ceil(0x10000) == 16, "ld_ceil 65536 = 16");
  BOOST_TEST(trng::int_math::log2_ceil(0x10001) == 17, "ld_ceil 65537 = 17");
}

BOOST_AUTO_TEST_CASE(test_ceil2) {
  BOOST_TEST(trng::int_math::ceil2(1) == 1, "ceil2 1 = 1");
  BOOST_TEST(trng::int_math::ceil2(0x10000) == 0x10000, "ceil 65536 = 65536");
  BOOST_TEST(trng::int_math::ceil2(0x10001) == 0x20000, "ceil 65537 = 131072");
}

BOOST_AUTO_TEST_CASE(test_mask) {
  BOOST_TEST(trng::int_math::mask(0x10000) == 0x1ffff, "mask 65536 = 65536");
  BOOST_TEST(trng::int_math::mask(0x10001) == 0x1ffff, "mask 65537 = 131071");
  BOOST_TEST(trng::int_math::mask(trng::int32_t(0x7fffffff)) == 0x7fffffff,
             "mask 2147483647 = 2147483647");
  BOOST_TEST(trng::int_math::mask(trng::uint32_t(0xffffffffu)) == 0xffffffffu,
             "mask 4294967295 = 4294967295");
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_linear_algebra)

BOOST_AUTO_TEST_CASE(test_matrix_vec_mult) {
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
  bool ok{true};
  for (int i{0}; i < n; ++i)
    ok = ok and (c[i] == c_exact[i]);
  BOOST_TEST(ok, "modular matrix-vector multiplication");
}

BOOST_AUTO_TEST_CASE(test_matrix_mult) {
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
  bool ok{true};
  for (int i{0}; i < n * n; ++i)
    ok = ok and (C[i] == C_exact[i]);
  BOOST_TEST(ok, "modular matrix-matrix multiplication");
}

BOOST_AUTO_TEST_CASE(test_gauss) {
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
  bool ok{true};
  for (int i{0}; i < n; ++i)
    ok = ok and (b[i] == x_exact[i]);
  BOOST_TEST(ok, "modular Gaussian elimination");
}

BOOST_AUTO_TEST_CASE(test_gauss_singular_matrix) {
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
  bool ok{true};
  for (int i{0}; i < n; ++i)
    ok = ok and (c[i] == b_orig[i]);
  BOOST_TEST(ok, "modular Gaussian elimination");
}

BOOST_AUTO_TEST_CASE(test_gauss_singular_matrix_no_solution) {
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
  BOOST_CHECK_THROW(trng::int_math::gauss<n>(A, b, m), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_modulo_inverse)

BOOST_AUTO_TEST_CASE(test_modulo_inverse_prime_modulus) {
  const long m{104729};  // must be prime
  bool modulo_inverse_ok{true};
  for (long a{1}; a < m and modulo_inverse_ok; ++a) {
    int b{trng::int_math::modulo_inverse(a, m)};
    modulo_inverse_ok = (static_cast<long long>(a) * static_cast<long long>(b)) % m == 1;
  }
  BOOST_TEST(modulo_inverse_ok, "modulo inverse prime modulus");
  BOOST_CHECK_THROW(trng::int_math::modulo_inverse(0, m), std::exception);
}

BOOST_AUTO_TEST_CASE(test_modulo_inverse_power2_modulus) {
  const long m{1024 * 1024};  // must be power of 2
  bool modulo_inverse_ok{true};
  for (long a{1}; a < m and modulo_inverse_ok; a += 2) {
    int b{trng::int_math::modulo_inverse(a, m)};
    modulo_inverse_ok = (static_cast<long long>(a) * static_cast<long long>(b)) % m == 1;
  }
  BOOST_TEST(modulo_inverse_ok, "modulo inverse power-of-2 modulus");
  BOOST_CHECK_THROW(trng::int_math::modulo_inverse(100, m), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
