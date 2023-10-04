// Copyright (c) 2000-2022, Heiko Bauke
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

#include <catch2/catch.hpp>

#include <trng/linear_algebra.hpp>

TEST_CASE("basic vector operations") {
  SECTION("comparison operations") {
    const trng::vector<int, 6> v1{1, 2, 3, 4, 5, 6};
    const trng::vector<int, 6> v2{1, 2, 3, 4, 5, 6};
    const trng::vector<int, 6> v3{1, 2, 3, 4, 5, 255};
    REQUIRE(v1 == v2);
    REQUIRE(v1 != v3);
  }
  SECTION("vector sum") {
    auto init1 = [](size_t i) { return std::uint8_t(i); };
    auto init2 = [](size_t i) { return std::uint8_t(3 * i); };
    const trng::vector<std::uint8_t, 256> v_256_1(init1);
    const trng::vector<std::uint8_t, 256> v_256_2(init2);
    const trng::vector<std::uint8_t, 256> v_256_3(v_256_1 + v_256_2);
    for (std::size_t i{0}; i < 256; ++i)
      REQUIRE(static_cast<std::uint8_t>(v_256_1(i) + v_256_2(i)) == v_256_3(i));
  }
  SECTION("scalar vector product") {
    auto init = [](size_t i) { return std::uint8_t(i); };
    const trng::vector<std::uint8_t, 256> v_256_1(init);
    const trng::vector<std::uint8_t, 256> v_256_2(std::uint8_t(13) * v_256_1);
    for (std::size_t i{0}; i < 256; ++i)
      REQUIRE(static_cast<std::uint8_t>(13 * v_256_1(i)) == v_256_2(i));
  }
  SECTION("vector scalar product") {
    auto init = [](size_t i) { return std::uint8_t(i); };
    const trng::vector<std::uint8_t, 256> v_256_1(init);
    const trng::vector<std::uint8_t, 256> v_256_2(v_256_1 * std::uint8_t(13));
    for (std::size_t i{0}; i < 256; ++i) {
      REQUIRE(static_cast<std::uint8_t>(13 * v_256_1(i)) == v_256_2(i));
    }
  }
}


TEST_CASE("basic matrix operations") {
  SECTION("matrix vector product") {
    const trng::matrix<int, 2> A{1, 3,  //
                                 2, 4};
    const trng::vector<int, 2> b{3, 7};
    const trng::vector<int, 2> c{24, 34};
    const auto c2{A * b};
    REQUIRE(c2 == c);
  }
  SECTION("matrix matrix product") {
    const trng::matrix<int, 2> A{1, 3,  //
                                 2, 4};
    const trng::matrix<int, 2> B{1, 3,  //
                                 2, -4};
    const trng::matrix<int, 2> C{7, -9,  //
                                 10, -10};
    const auto C2{A * B};
    REQUIRE(C2 == C);
  }
  SECTION("matrix power") {
    const trng::matrix<int, 2> A{1, 3,  //
                                 2, 4};
    const trng::matrix<int, 2> A_5{1069, 2337,  //
                                   1558, 3406};
    const auto A_5_2{trng::power(A, 5)};
    REQUIRE(A_5_2 == A_5);
  }
}


TEST_CASE("GF(2)") {
  SECTION("addition") {
    const trng::GF2 zero(false);
    const trng::GF2 one(true);
    REQUIRE(zero + zero == zero);
    REQUIRE(zero + one == one);
    REQUIRE(one + zero == one);
    REQUIRE(one + one == zero);
  }

  SECTION("multiplication") {
    const trng::GF2 zero(false);
    const trng::GF2 one(true);
    REQUIRE(zero * zero == zero);
    REQUIRE(zero * one == zero);
    REQUIRE(one * zero == zero);
    REQUIRE(one * one == one);
  }

  SECTION("matrix power") {
    const trng::GF2 zero(false);
    const trng::GF2 one(true);
    const trng::matrix<trng::GF2, 4> A{one,  one, zero, one,   //
                                       one,  one, one,  zero,  //
                                       one,  one, zero, zero,  //
                                       zero, one, zero, one};
    const trng::matrix<trng::GF2, 4> A_8{one,  one, zero, one,   //
                                         one,  one, one,  zero,  //
                                         one,  one, zero, zero,  //
                                         zero, one, zero, one};
    auto A_8_2{trng::power(A, 8)};
    REQUIRE(A_8_2 == A_8);
  }
}
