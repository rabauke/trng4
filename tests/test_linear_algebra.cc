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

#include <cinttypes>
#include <ciso646>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/linear_algebra.hpp>


//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_linear_algebra)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_vector)
BOOST_AUTO_TEST_CASE(test_basic) {
  {
    trng::vector<int, 6> v1{1, 2, 3, 4, 5, 6};
    trng::vector<int, 6> v2{1, 2, 3, 4, 5, 255};
    {
      const bool ok = v1 == v1;
      BOOST_TEST(ok, "vectors are equal");
    }
    {
      const bool ok = v1 != v2;
      BOOST_TEST(ok, "vectors do not equal");
    }
  }
  {
    auto init1 = [](size_t i) { return std::uint8_t(i); };
    auto init2 = [](size_t i) { return std::uint8_t(3 * i); };
    trng::vector<std::uint8_t, 256> v_256_1(init1);
    trng::vector<std::uint8_t, 256> v_256_2(init2);
    trng::vector<std::uint8_t, 256> v_256_3(v_256_1 + v_256_2);
    bool ok{true};
    for (std::size_t i{0}; i < 256; ++i) {
      ok = ok and (static_cast<std::uint8_t>(v_256_1(i) + v_256_2(i)) == v_256_3(i));
    }
    BOOST_TEST(ok, "sum vector ok");
  }
  {
    auto init = [](size_t i) { return std::uint8_t(i); };
    trng::vector<std::uint8_t, 256> v_256_1(init);
    trng::vector<std::uint8_t, 256> v_256_2(std::uint8_t(13) * v_256_1);
    bool ok{true};
    for (std::size_t i{0}; i < 256; ++i) {
      ok = ok and (static_cast<std::uint8_t>(13 * v_256_1(i)) == v_256_2(i));
    }
    BOOST_TEST(ok, "product vector ok");
  }
  {
    auto init = [](size_t i) { return std::uint8_t(i); };
    trng::vector<std::uint8_t, 256> v_256_1(init);
    trng::vector<std::uint8_t, 256> v_256_2(v_256_1 * std::uint8_t(13));
    bool ok{true};
    for (std::size_t i{0}; i < 256; ++i) {
      ok = ok and (static_cast<std::uint8_t>(13 * v_256_1(i)) == v_256_2(i));
    }
    BOOST_TEST(ok, "product vector ok");
  }
  {
    trng::matrix<int, 2> A{1, 3,  //
                           2, 4};
    trng::vector<int, 2> b{3, 7};
    trng::vector<int, 2> c{24, 34};
    auto c2 = A * b;
    const bool ok{c2 == c};
    BOOST_TEST(ok, "matrix vector product ok");
  }
  {
    trng::matrix<int, 2> A{1, 3,  //
                           2, 4};
    trng::matrix<int, 2> B{1, 3,  //
                           2, -4};
    trng::matrix<int, 2> C{7, -9,  //
                           10, -10};
    auto C2 = A * B;
    const bool ok{C2 == C};
    BOOST_TEST(ok, "matrix matrix product ok");
  }
  {
    trng::matrix<int, 2> A{1, 3,  //
                           2, 4};
    trng::matrix<int, 2> A_5{1069, 2337,  //
                             1558, 3406};
    auto A_5_2 = trng::power(A, 5);
    const bool ok{A_5_2 == A_5};
    BOOST_TEST(ok, "matrix matrix power ok");
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_GF2)
BOOST_AUTO_TEST_CASE(test_add) {
  const trng::GF2 zero(false);
  const trng::GF2 one(true);
  BOOST_TEST(zero + zero == zero, "addition 0 + 0");
  BOOST_TEST(zero + one == one, "addition 0 + 1");
  BOOST_TEST(one + zero == one, "addition 1 + 0");
  BOOST_TEST(one + one == zero, "addition 1 + 1");
}

BOOST_AUTO_TEST_CASE(test_mult) {
  const trng::GF2 zero(false);
  const trng::GF2 one(true);
  BOOST_TEST((zero * zero == zero), "multiplication 0 * 0");
  BOOST_TEST((zero * one == zero), "multiplication 0 * 1");
  BOOST_TEST((one * zero == zero), "multiplication 1 * 0");
  BOOST_TEST((one * one == one), "multiplication 1 * 1");
}

BOOST_AUTO_TEST_CASE(test_matrix_power) {
  const trng::GF2 zero(false);
  const trng::GF2 one(true);
  trng::matrix<trng::GF2, 4> A{one,  one, zero, one,   //
                               one,  one, one,  zero,  //
                               one,  one, zero, zero,  //
                               zero, one, zero, one};
  trng::matrix<trng::GF2, 4> A_8{one,  one, zero, one,   //
                                 one,  one, one,  zero,  //
                                 one,  one, zero, zero,  //
                                 zero, one, zero, one};
  auto A_8_2 = trng::power(A, 8);
  const bool ok{A_8_2 == A_8};
  BOOST_TEST(ok, "matrix matrix power in GF2 ok");
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
