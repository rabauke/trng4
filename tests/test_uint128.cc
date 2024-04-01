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


#include <catch2/catch.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#if defined(__SIZEOF_INT128__)
#undef __SIZEOF_INT128__
#endif
#include <trng/uint128.hpp>

template<typename TA = trng::uint128, typename TY = trng::uint128>
struct unary_op_tuple {
  TA a;
  TY y;
};


template<typename TA = trng::uint128, typename TB = trng::uint128, typename TY = trng::uint128>
struct binary_op_tuple {
  TA a;
  TB b;
  TY y;
};


template<typename TA = trng::uint128, typename TB = trng::uint128>
struct pair {
  TA a;
  TB b;
};


TEST_CASE("uint128 cast") {
  SECTION("to uint64") {
    using data_type = unary_op_tuple<trng::uint128, std::uint64_t>;
    // clang-format off
    auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, std::uint64_t{0x0}},
        data_type{trng::uint128{0x2, 0x1}, std::uint64_t{0x1}},
        data_type{trng::uint128{0x3, 0xffffffff}, std::uint64_t{0xffffffff}},
        data_type{trng::uint128{0x4, 0xffffffffffffffff}, std::uint64_t{0xffffffffffffffff}}
    )};
    // clang-format on
    trng::uint128 a{a_y.a};
    std::uint64_t y{a_y.y};
    REQUIRE(static_cast<std::uint64_t>(a) == y);
  }
  SECTION("to float") {
    using data_type = unary_op_tuple<trng::uint128, float>;
    // clang-format off
    auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, 0.f},
        data_type{trng::uint128{0x2, 0x1}, 36893488147419103233.f},
        data_type{trng::uint128{0x3, 0xffffffffffffffff}, 73786976294838206463.f}
        // 340282366920938463463374607431768211455 is larger than FLT_MAX on most systems
        // data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, 340282366920938463463374607431768211455.f}
    )};
    // clang-format on
    trng::uint128 a{a_y.a};
    float y{a_y.y};
    REQUIRE(static_cast<float>(a) == y);
  }
  SECTION("to double") {
    using data_type = unary_op_tuple<trng::uint128, double>;
    // clang-format off
    auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, 0.},
        data_type{trng::uint128{0x2, 0x1}, 36893488147419103233.},
        data_type{trng::uint128{0x3, 0xffffffffffffffff}, 73786976294838206463.},
        data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, 340282366920938463463374607431768211455.}
    )};
    // clang-format on
    trng::uint128 a{a_y.a};
    double y{a_y.y};
    REQUIRE(static_cast<double>(a) == y);
  }
  SECTION("to long double") {
    using data_type = unary_op_tuple<trng::uint128, long double>;
    // clang-format off
    auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, 0.},
        data_type{trng::uint128{0x2, 0x1}, 36893488147419103233.l},
        data_type{trng::uint128{0x3, 0xffffffffffffffff}, 73786976294838206463.l},
        data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, 340282366920938463463374607431768211455.l}
    )};
    // clang-format on
    trng::uint128 a{a_y.a};
    long double y{a_y.y};
    REQUIRE(static_cast<long double>(a) == y);
  }
}


TEST_CASE("uint128 math") {
  SECTION("arithmetic operators") {
    SECTION("unary plus") {
      using data_type = unary_op_tuple<>;
      // clang-format off
      auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, trng::uint128{0x0, 0x0}},
        data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}},
        data_type{trng::uint128{0xa4f2087e8b36d663, 0x60ec8b7d23e07770}, trng::uint128{0xa4f2087e8b36d663, 0x60ec8b7d23e07770}},
        data_type{trng::uint128{0xcfc3bcb6332ade50, 0x886741fcf2ecc125}, trng::uint128{0xcfc3bcb6332ade50, 0x886741fcf2ecc125}},
        data_type{trng::uint128{0x59c94057fd77d893, 0xf062680f3d759e10}, trng::uint128{0x59c94057fd77d893, 0xf062680f3d759e10}},
        data_type{trng::uint128{0xe03dc9c3043ab169, 0x8cd70bce637c48af}, trng::uint128{0xe03dc9c3043ab169, 0x8cd70bce637c48af}},
        data_type{trng::uint128{0x9f5a86f9c799b8ed, 0x96488fc27a04f3e1}, trng::uint128{0x9f5a86f9c799b8ed, 0x96488fc27a04f3e1}}
      )};
      // clang-format on
      trng::uint128 a{a_y.a};
      trng::uint128 y{a_y.y};
      REQUIRE(+a == y);
    }
    SECTION("unary minus") {
      using data_type = unary_op_tuple<>;
      // clang-format off
      auto a_y{GENERATE(
        data_type{trng::uint128{0x0, 0x0}, trng::uint128{0x0, 0x0}},
        data_type{trng::uint128{0x0, 0x1}, trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}},
        data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, trng::uint128{0x0, 0x1}},
        data_type{trng::uint128{0x352e99558a467bfe, 0x4d1e7adb04b931c9}, trng::uint128{0xcad166aa75b98401, 0xb2e18524fb46ce37}},
        data_type{trng::uint128{0x788e6ea0e67267f5, 0x3b0145b631ca7310}, trng::uint128{0x8771915f198d980a, 0xc4feba49ce358cf0}},
        data_type{trng::uint128{0xd8249eece2d0c5d1, 0x31dc129870f4055d}, trng::uint128{0x27db61131d2f3a2e, 0xce23ed678f0bfaa3}},
        data_type{trng::uint128{0x39d6029725d81f3b, 0x2bb385c93d500890}, trng::uint128{0xc629fd68da27e0c4, 0xd44c7a36c2aff770}},
        data_type{trng::uint128{0x3a48350dccb6120b, 0x6f89f2f9791783d6}, trng::uint128{0xc5b7caf23349edf4, 0x90760d0686e87c2a}}
      )};
      // clang-format on
      trng::uint128 a{a_y.a};
      trng::uint128 y{a_y.y};
      REQUIRE(-a == y);
    }
    SECTION("plus") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0xd091bb5c22ae9ef6, 0xe7e1faeed5c31f79}, trng::uint128{0xa7de9f4ccc450cba, 0x924668f5c7dc380}, trng::uint128{0x78705aa8eef3abb0, 0xf106617e3240e2f9}},
        data_type{trng::uint128{0x2082352cf807b7df, 0xe9d300053895afe1}, trng::uint128{0xd96089c53640ac4c, 0xef1a2e6dae6d9426}, trng::uint128{0xf9e2bef22e48642c, 0xd8ed2e72e7034407}},
        data_type{trng::uint128{0xa1e24bba4ee4092b, 0x18f868638c16a625}, trng::uint128{0xadc1965b6613ba46, 0xc1fb41c2bd9b0ecd}, trng::uint128{0x4fa3e215b4f7c371, 0xdaf3aa2649b1b4f2}},
        data_type{trng::uint128{0x474ba8c43039cd1a, 0x8c006d5ffe2d7810}, trng::uint128{0xbe3dedfc7989c8ee, 0x6468fd6e6c0df032}, trng::uint128{0x58996c0a9c39608, 0xf0696ace6a3b6842}},
        data_type{trng::uint128{0xf51f2ae7ff1816e4, 0xf702ef59f7badafa}, trng::uint128{0xa7cd66342c826d8b, 0x2bd2e4124d4a2dbe}, trng::uint128{0x9cec911c2b9a8470, 0x22d5d36c450508b8}},
        data_type{trng::uint128{0x285954a1b9d09511, 0xf878c4b3fb2a0137}, trng::uint128{0xb4bf6fa7cc1a8959, 0x826328251097330}, trng::uint128{0xdd18c44985eb1e6b, 0x9ef7364c337467}},
        data_type{trng::uint128{0xf508e4aa1c1fe652, 0x7c419418cc50aa59}, trng::uint128{0x46e46cb0df577ec2, 0xbd1e364262c5564}, trng::uint128{0x3bed515afb776514, 0x8813777cf27cffbd}},
        data_type{trng::uint128{0xccdf2e5c4c0a1f3b, 0x2452a9dc01397d8d}, trng::uint128{0x18dda0c9fe7b45d9, 0xd2ce21c9d268409a}, trng::uint128{0xe5bccf264a856514, 0xf720cba5d3a1be27}},
        data_type{trng::uint128{0x6bf88c311cca797a, 0xea6da4aea3c78807}, trng::uint128{0xb1e049e1200bfa47, 0x512d6e73c3851eee}, trng::uint128{0x1dd8d6123cd673c2, 0x3b9b1322674ca6f5}},
        data_type{trng::uint128{0xcace1969e0e0d4ad, 0xf5a14bab80f00988}, trng::uint128{0xf341c0817d973e48, 0x8d17554a9e20d28}, trng::uint128{0xbe0fd9eb5e7812f5, 0xfe72c1002ad216b0}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a + b == y);
      a += b;
      REQUIRE(a == y);
    }
    SECTION("minus") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0x70518ce6203ac303, 0x61add0ab35d0430c}, trng::uint128{0xc05309bed23d2d63, 0x414de9c5d2229f23}, trng::uint128{0xaffe83274dfd95a0, 0x205fe6e563ada3e9}},
        data_type{trng::uint128{0xc3f8e8920d1c8509, 0xcb92388e095436bf}, trng::uint128{0x818666a3f0a8b109, 0xb2f6b12769a48341}, trng::uint128{0x427281ee1c73d400, 0x189b87669fafb37e}},
        data_type{trng::uint128{0x2fd6e20868a29af9, 0x7d61330b753ec6fc}, trng::uint128{0xe4123c566c548c8f, 0xf5941f6194b993aa}, trng::uint128{0x4bc4a5b1fc4e0e69, 0x87cd13a9e0853352}},
        data_type{trng::uint128{0x7211efea7cd15133, 0xa574c4ffcb41f198}, trng::uint128{0x8c1651342876763c, 0x237ce42ec300d11b}, trng::uint128{0xe5fb9eb6545adaf7, 0x81f7e0d10841207d}},
        data_type{trng::uint128{0xb598eef6ebbe7347, 0xc1332568ceba5a70}, trng::uint128{0x263821ca3aeb8202, 0x41ec0f84cf4ac36d}, trng::uint128{0x8f60cd2cb0d2f145, 0x7f4715e3ff6f9703}},
        data_type{trng::uint128{0x46a99459b4ad9f11, 0xae00feaa00b8b573}, trng::uint128{0xd7393ee6fd0fc06a, 0x4118a30a551b54a4}, trng::uint128{0x6f705572b79ddea7, 0x6ce85b9fab9d60cf}},
        data_type{trng::uint128{0xa7b480b6b5f0b06c, 0x29a0ec27a4daa010}, trng::uint128{0xd074f86f4cc1c54a, 0x3e57a70303774cda}, trng::uint128{0xd73f8847692eeb21, 0xeb494524a1635336}},
        data_type{trng::uint128{0x1e76a1c574be9133, 0x7f94c950c61f6ed6}, trng::uint128{0xede43895379ce627, 0x59988939e8490ddc}, trng::uint128{0x309269303d21ab0c, 0x25fc4016ddd660fa}},
        data_type{trng::uint128{0xf5b1c7a192e195f8, 0x572384d4e0732c88}, trng::uint128{0x325410e1d9352f6a, 0x4047080af47c081d}, trng::uint128{0xc35db6bfb9ac668e, 0x16dc7cc9ebf7246b}},
        data_type{trng::uint128{0x95d41b68cee496c3, 0x394bbd52048cd47c}, trng::uint128{0x9db51a85c765d71f, 0x79297527fcca2773}, trng::uint128{0xf81f00e3077ebfa3, 0xc022482a07c2ad09}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a - b == y);
      a -= b;
      REQUIRE(a == y);
    }
    SECTION("minus") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0x70518ce6203ac303, 0x61add0ab35d0430c}, trng::uint128{0xc05309bed23d2d63, 0x414de9c5d2229f23}, trng::uint128{0xaffe83274dfd95a0, 0x205fe6e563ada3e9}},
        data_type{trng::uint128{0xc3f8e8920d1c8509, 0xcb92388e095436bf}, trng::uint128{0x818666a3f0a8b109, 0xb2f6b12769a48341}, trng::uint128{0x427281ee1c73d400, 0x189b87669fafb37e}},
        data_type{trng::uint128{0x2fd6e20868a29af9, 0x7d61330b753ec6fc}, trng::uint128{0xe4123c566c548c8f, 0xf5941f6194b993aa}, trng::uint128{0x4bc4a5b1fc4e0e69, 0x87cd13a9e0853352}},
        data_type{trng::uint128{0x7211efea7cd15133, 0xa574c4ffcb41f198}, trng::uint128{0x8c1651342876763c, 0x237ce42ec300d11b}, trng::uint128{0xe5fb9eb6545adaf7, 0x81f7e0d10841207d}},
        data_type{trng::uint128{0xb598eef6ebbe7347, 0xc1332568ceba5a70}, trng::uint128{0x263821ca3aeb8202, 0x41ec0f84cf4ac36d}, trng::uint128{0x8f60cd2cb0d2f145, 0x7f4715e3ff6f9703}},
        data_type{trng::uint128{0x46a99459b4ad9f11, 0xae00feaa00b8b573}, trng::uint128{0xd7393ee6fd0fc06a, 0x4118a30a551b54a4}, trng::uint128{0x6f705572b79ddea7, 0x6ce85b9fab9d60cf}},
        data_type{trng::uint128{0xa7b480b6b5f0b06c, 0x29a0ec27a4daa010}, trng::uint128{0xd074f86f4cc1c54a, 0x3e57a70303774cda}, trng::uint128{0xd73f8847692eeb21, 0xeb494524a1635336}},
        data_type{trng::uint128{0x1e76a1c574be9133, 0x7f94c950c61f6ed6}, trng::uint128{0xede43895379ce627, 0x59988939e8490ddc}, trng::uint128{0x309269303d21ab0c, 0x25fc4016ddd660fa}},
        data_type{trng::uint128{0xf5b1c7a192e195f8, 0x572384d4e0732c88}, trng::uint128{0x325410e1d9352f6a, 0x4047080af47c081d}, trng::uint128{0xc35db6bfb9ac668e, 0x16dc7cc9ebf7246b}},
        data_type{trng::uint128{0x95d41b68cee496c3, 0x394bbd52048cd47c}, trng::uint128{0x9db51a85c765d71f, 0x79297527fcca2773}, trng::uint128{0xf81f00e3077ebfa3, 0xc022482a07c2ad09}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a - b == y);
      a -= b;
      REQUIRE(a == y);
    }
    SECTION("multiply") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0x5a065b97114dee4f, 0xd4b12f5fcb29360a}, trng::uint128{0x2984c787ed702bbe, 0xcb563b4d6fa56696}, trng::uint128{0x24763a4b1d26603, 0xfabc4c13a01fa5dc}},
        data_type{trng::uint128{0x95d3de16983162a8, 0x8cbaafb3bb98b27f}, trng::uint128{0x4fabc9eddcd87a48, 0x874df2959ecfe9f0}, trng::uint128{0xe2536339b924c138, 0xadb2ecb904dee10}},
        data_type{trng::uint128{0xeacd3439b1fac842, 0x492cbef1ae08ab78}, trng::uint128{0x2a67f49f1e9aa4e1, 0x9a1b7d0878d22934}, trng::uint128{0xa6f5e3f20c286b74, 0xe5274533a5a90c60}},
        data_type{trng::uint128{0xc1d7dfd0646f1d40, 0xc0f463c48fc23a81}, trng::uint128{0x435216025718a361, 0xa771ba4487a3b97c}, trng::uint128{0xc73f32843d47a76, 0x887a98e2457e8f7c}},
        data_type{trng::uint128{0x6164e6233543f2bc, 0x915cc2538701d0df}, trng::uint128{0xb0705c82b7526048, 0xbf86dcd7fd066ea4}, trng::uint128{0x5696e5022a86fc29, 0xbd0e66458d23a0dc}},
        data_type{trng::uint128{0x136b2fdd677a359e, 0xdcfacd05a4ea31e}, trng::uint128{0x7356b1bbb872426d, 0x1575515de99eadb3}, trng::uint128{0x89d1f2029c7474aa, 0xfdabb19b43bb53fa}},
        data_type{trng::uint128{0x87e2593597c34e42, 0xc77780f05b396fba}, trng::uint128{0x3a9e3c0f8168599c, 0xe9d07a328eeab382}, trng::uint128{0x1f580c875017eb41, 0x49bb3ea4c84dca74}},
        data_type{trng::uint128{0xef1b52e6f7080941, 0x2141888b278946b0}, trng::uint128{0x27023ee880d10fac, 0xd368bdc27664b5a7}, trng::uint128{0x659fb28888348175, 0x9896af4f96478cd0}},
        data_type{trng::uint128{0x919e6d646518b459, 0x7829fc226325d30e}, trng::uint128{0x89d0cf468bed7368, 0xff02af497294e430}, trng::uint128{0x947f48a83fc8d24, 0xb9a1d19887280aa0}},
        data_type{trng::uint128{0x30c0399ba19b463, 0x564dab7563794f97}, trng::uint128{0x14034fbbdabd4cc4, 0x71535cf89aaeea20}, trng::uint128{0xa6f23c77915994c7, 0x5b23b036408bf8e0}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a * b == y);
      a *= b;
      REQUIRE(a == y);
    }
    SECTION("divide") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0x1b4d989d7fa09780, 0xf63ef3d2fadc6788}, trng::uint128{0xda603f4888a, 0xfd7149f3f014d704}, trng::uint128{0x0, 0x2001d}},
        data_type{trng::uint128{0x12fb56808c904fa, 0xc660883ffa1cce2a}, trng::uint128{0x59d803fb331, 0x5c9dc83645273235}, trng::uint128{0x0, 0x3616}},
        data_type{trng::uint128{0xd13ac8b85cf9c9b3, 0xde62c6bdadf500ad}, trng::uint128{0x66dce43761e, 0x924413728d90c32}, trng::uint128{0x0, 0x208b86}},
        data_type{trng::uint128{0x159d967e58a2c06c, 0x665827cbdb1aa208}, trng::uint128{0x3d6b25290a5, 0xc50941f91d46440c}, trng::uint128{0x0, 0x5a18a}},
        data_type{trng::uint128{0x4286ddf10b8905b4, 0xccd149a4a8fd9757}, trng::uint128{0x2f1494b9813, 0x8ac023d6d8755d29}, trng::uint128{0x0, 0x169bd7}},
        data_type{trng::uint128{0x6e7122f0bffc21b1, 0xe9203368220c0724}, trng::uint128{0x6ad203b6fb5, 0x234e00cb62703d2c}, trng::uint128{0x0, 0x108ade}},
        data_type{trng::uint128{0x2e8d86cbfb7bfb5e, 0x438896871869325b}, trng::uint128{0xe718672d4d4, 0x48df1f1dd92d70c4}, trng::uint128{0x0, 0x3391d}},
        data_type{trng::uint128{0x25420afc485d46db, 0x22d56381cd572d60}, trng::uint128{0x7da944f13f7, 0x2f0c37d3fa9308b4}, trng::uint128{0x0, 0x4be72}},
        data_type{trng::uint128{0xde89ef2b13dac708, 0x9467851da09c428d}, trng::uint128{0x5674c771e34, 0x20514e669edd5580}, trng::uint128{0x0, 0x292f22}},
        data_type{trng::uint128{0x8cc3a36c0212714e, 0x251bc1f6ae274af0}, trng::uint128{0x5e86f504088, 0x449f41c77c8a029b}, trng::uint128{0x0, 0x17d384}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a / b == y);
      a /= b;
      REQUIRE(a == y);
    }
    SECTION("reminder") {
      using data_type = binary_op_tuple<>;
      // clang-format off
      auto a_b_y{GENERATE(
        data_type{trng::uint128{0xc7be9961e09aebe7, 0x63c5ecb935d657e1}, trng::uint128{0x8c08c64db7e, 0xda5894bdbae3349a}, trng::uint128{0xc2cadb34c7, 0xa7165e947c867b45}},
        data_type{trng::uint128{0x3ddf7ae445d3249c, 0x6766c9407e10ad9c}, trng::uint128{0xbea4300c9f9, 0x97956305f5b0ccfd}, trng::uint128{0xf2d6010901, 0xcf0fdf9b64a8c8b3}},
        data_type{trng::uint128{0x18b13e7f39320ca2, 0x21c900787d84661c}, trng::uint128{0xafd083db500, 0x85e02efa964462b}, trng::uint128{0x8c377809ad4, 0xc895ee46eeab73db}},
        data_type{trng::uint128{0xf12a3a21f4772b41, 0xf4c53baba6e76b3b}, trng::uint128{0x5e5501a5580, 0xbc7bba02889add0d}, trng::uint128{0x3d58453c760, 0x80447095dc3313e7}},
        data_type{trng::uint128{0x9340ded2c1ebc21f, 0xf4db6546f6c42e3}, trng::uint128{0xc7bcfc3cab5, 0x6454f14c4a882612}, trng::uint128{0x916847bdf14, 0x2e4a4b08b9435e4d}},
        data_type{trng::uint128{0x3c1a89438d899f74, 0x5a6899ba0d9b6827}, trng::uint128{0xedec3d2f756, 0xd88c6951b284db}, trng::uint128{0x755ef1382e3, 0x16c9bda2f6fcd7e4}},
        data_type{trng::uint128{0xd239c5cb5290106a, 0x3f17adb67acdc26}, trng::uint128{0x7c9e604de91, 0xac9746f946da27d1}, trng::uint128{0x1c37e22c5a8, 0x2f6f6bbf58a14a95}},
        data_type{trng::uint128{0xb039b90e88f1afa, 0x2b42ee31cf239e4d}, trng::uint128{0x72606b78c17, 0xb72574e6ceb4e1ba}, trng::uint128{0x2e3881468e9, 0xd33f78e29a025c61}},
        data_type{trng::uint128{0xa62c6e93421fae11, 0xbb522891213d9f32}, trng::uint128{0x822d9f9a64c, 0x24df682c1eccbdec}, trng::uint128{0x719c7aa2f56, 0x2ba4523e3b1d00ae}},
        data_type{trng::uint128{0xa5d2adeb7e4aab21, 0x736fbc7560e56773}, trng::uint128{0xd1500cb65c6, 0x15e22cb7a1247eac}, trng::uint128{0x2c04df82740, 0x735402bdb3cc9cd7}}
      )};
      // clang-format on
      trng::uint128 a{a_b_y.a};
      trng::uint128 b{a_b_y.b};
      trng::uint128 y{a_b_y.y};
      REQUIRE(a % b == y);
      a %= b;
      REQUIRE(a == y);
    }
  }
}


TEST_CASE("uint128 comparison operators") {
  SECTION("less than") {
    using data_type = pair<>;
    // clang-format off
    auto a_b{GENERATE(
      data_type{trng::uint128{0x4fa4664e16fea67f, 0xec629bbfa5014386}, trng::uint128{0x984b1bef73161b54, 0x43204f200ac40f25}},
      data_type{trng::uint128{0x6e221244a21075d9, 0x2f501f52959a12f4}, trng::uint128{0x9a53eca9cc200dd8, 0xb6123cc02ae4efad}},
      data_type{trng::uint128{0x38c464c2d4ca75a4, 0x1e0f15595330cf6e}, trng::uint128{0xe7a64774ba060582, 0xfad0ca2f5ac9908f}},
      data_type{trng::uint128{0x4bf2f32ba7e130fd, 0x519b74626b919194}, trng::uint128{0x7059c853ae2f213e, 0x1c724f28b51305fc}},
      data_type{trng::uint128{0x42108734298e5c9d, 0x68a1dd2b223c8c36}, trng::uint128{0x6c963148c82b32bb, 0x820300241cfa2fa0}}
    )};
    // clang-format on
    trng::uint128 a{a_b.a};
    trng::uint128 b{a_b.b};
    REQUIRE(a < b);
  }
  SECTION("less equal than") {
    using data_type = pair<>;
    // clang-format off
    auto a_b{GENERATE(
      data_type{trng::uint128{0x15e45ce2fe584a91, 0x4332093f2e7b9117}, trng::uint128{0xf68cf8789b1e1c75, 0x8bfb75ae39ff446d}},
      data_type{trng::uint128{0x85692875309da59f, 0x3b49c50966cd636e}, trng::uint128{0xcd0f4b4dc34f792c, 0x77afaf544136041}},
      data_type{trng::uint128{0x7d28670869cd6a2e, 0x9fc266e32b8f1988}, trng::uint128{0xedc962973e7fd864, 0xbaf6f6ba19c9ff6e}},
      data_type{trng::uint128{0x7d157a585825dfde, 0x941a37e04818babb}, trng::uint128{0xaddbd3af4c8ea8e7, 0x65407c37a2ef9aa6}},
      data_type{trng::uint128{0x3cbe9dc1f7f8d0ce, 0x75771de936b9cf69}, trng::uint128{0x5e10541b26cd9065, 0xfcec6367c4ed1ef8}}
    )};
    // clang-format on
    trng::uint128 a{a_b.a};
    trng::uint128 b{a_b.b};
    REQUIRE(a <= b);
    REQUIRE(a <= a);
    REQUIRE(b <= b);
  }
  SECTION("greater than") {
    using data_type = pair<>;
    // clang-format off
    auto a_b{GENERATE(
      data_type{trng::uint128{0x1b54beda808164b1, 0xa75ca4570068b861}, trng::uint128{0x9a9410fe24bc427, 0xe29a5eddf58f8c10}},
      data_type{trng::uint128{0xe9cd2a63049ffdb7, 0xcbd2b4cc95356bb0}, trng::uint128{0x7e822ee10335be36, 0xc76fef0e495dbce8}},
      data_type{trng::uint128{0xb70cab9d19445725, 0xe75a3b16627f5136}, trng::uint128{0x19453535c4508ca2, 0x4309fd7e53ea8de9}},
      data_type{trng::uint128{0xe4137f415af821c1, 0x558bb5a6b85003e7}, trng::uint128{0x55d9f238210a7aea, 0xae02a6ab4abdf123}},
      data_type{trng::uint128{0xb2e101c210c101fc, 0x32a3aa3d838c4690}, trng::uint128{0x22f5256bd8dc2d8e, 0xb8a25d9d3b13600f}}
    )};
    // clang-format on
    trng::uint128 a{a_b.a};
    trng::uint128 b{a_b.b};
    REQUIRE(a > b);
  }
  SECTION("greater equal than") {
    using data_type = pair<>;
    // clang-format off
    auto a_b{GENERATE(
      data_type{trng::uint128{0x2ed401c83eb994dd, 0x3d6c2f3d553fd4da}, trng::uint128{0x7d188006c89d813, 0xbe7ba68db096d0b9}},
      data_type{trng::uint128{0xe2f27231809ad218, 0x757306e989740b3}, trng::uint128{0x8001786a41026a58, 0x7adc2d6566c9fc5f}},
      data_type{trng::uint128{0xe79be06843b4bc3b, 0x9c203868ef2b0be2}, trng::uint128{0x7d6a2d714f66622b, 0x2afd45f7f687016b}},
      data_type{trng::uint128{0xfa8ad0b33b7e96f2, 0xb6732508bf351d33}, trng::uint128{0x9e1f626778f76a0a, 0xdc04692d9c6500f3}},
      data_type{trng::uint128{0xce348e127d07e7ee, 0x93a4057983b78de9}, trng::uint128{0x801ee8985291851f, 0x78993f7e5177a95e}}
    )};
    // clang-format on
    trng::uint128 a{a_b.a};
    trng::uint128 b{a_b.b};
    REQUIRE(a >= b);
    REQUIRE(a >= a);
    REQUIRE(b >= b);
  }
}


TEST_CASE("uint128 bitwise shift operators") {
  SECTION("bitwise shift operators left") {
    using data_type = binary_op_tuple<trng::uint128, int, trng::uint128>;
    // clang-format off
    auto a_b_y{GENERATE(
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 127, trng::uint128{0x8000000000000000, 0x0}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 72, trng::uint128{0x95b5e088e37a7b00, 0x0}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 65, trng::uint128{0x5d2b6bc111c6f4f6, 0x0}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 64, trng::uint128{0xae95b5e088e37a7b, 0x0}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 63, trng::uint128{0x574adaf04471bd3d, 0x8000000000000000}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 8, trng::uint128{0x432e84f3d5350cae, 0x95b5e088e37a7b00}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, 0, trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -8, trng::uint128{0xf432e84f3d535, 0xcae95b5e088e37a}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -63, trng::uint128{0x0, 0x1e865d09e7aa6a19}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -64, trng::uint128{0x0, 0xf432e84f3d5350c}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -65, trng::uint128{0x0, 0x7a1974279ea9a86}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -72, trng::uint128{0x0, 0xf432e84f3d535}},
      data_type{trng::uint128{0xf432e84f3d5350c, 0xae95b5e088e37a7b}, -127, trng::uint128{0x0, 0x0}}
    )};
    // clang-format on
    trng::uint128 a{a_b_y.a};
    int b{a_b_y.b};
    trng::uint128 c{a_b_y.y};
    REQUIRE(a << b == c);
    a <<= b;
    REQUIRE(a == c);
  }
  SECTION("bitwise shift operators right") {
    using data_type = binary_op_tuple<trng::uint128, int, trng::uint128>;
    // clang-format off
    auto a_b_y{GENERATE(
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 127, trng::uint128{0x0, 0x0}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 72, trng::uint128{0x0, 0xadcc451df2e9f}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 65, trng::uint128{0x0, 0x56e6228ef974fad}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 64, trng::uint128{0x0, 0xadcc451df2e9f5b}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 63, trng::uint128{0x0, 0x15b988a3be5d3eb6}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 8, trng::uint128{0xadcc451df2e9f, 0x5b124a3fe8ef421e}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, 0, trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -8, trng::uint128{0xdcc451df2e9f5b12, 0x4a3fe8ef421e8800}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -63, trng::uint128{0x89251ff477a10f44, 0x0}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -64, trng::uint128{0x124a3fe8ef421e88, 0x0}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -65, trng::uint128{0x24947fd1de843d10, 0x0}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -72, trng::uint128{0x4a3fe8ef421e8800, 0x0}},
      data_type{trng::uint128{0xadcc451df2e9f5b, 0x124a3fe8ef421e88}, -127, trng::uint128{0x0, 0x0}}
    )};
    // clang-format on
    trng::uint128 a{a_b_y.a};
    int b{a_b_y.b};
    trng::uint128 c{a_b_y.y};
    REQUIRE(a >> b == c);
    a >>= b;
    REQUIRE(a == c);
  }
}


TEST_CASE("uint128 i/o operations") {
  SECTION("decimal output") {
    using data_type = unary_op_tuple<trng::uint128, std::string>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{trng::uint128{0x0, 0x0}, std::string {"0"}},
      data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, std::string{"340282366920938463463374607431768211455"}},
      data_type{trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}, std::string{"177508241696884716594379275141355476411"}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << std::dec << a_y.a;
      REQUIRE(strstr.str() == a_y.y);
    }
    {
      std::stringstream strstr;
      strstr << std::dec << std::showbase << a_y.a;
      REQUIRE(strstr.str() == a_y.y);
    }
  }
  SECTION("hex output") {
    using data_type = unary_op_tuple<trng::uint128, std::string>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{trng::uint128{0x0, 0x0}, std::string {"0"}},
      data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, std::string{"ffffffffffffffffffffffffffffffff"}},
      data_type{trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}, std::string{"858ad803aac7fd1118c34c954a2915bb"}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << std::hex << a_y.a;
      REQUIRE(strstr.str() == a_y.y);
    }
    {
      std::stringstream strstr;
      strstr << std::hex << std::showbase << a_y.a;
      REQUIRE(strstr.str() == "0x" + a_y.y);
    }
  }
  SECTION("octal output") {
    using data_type = unary_op_tuple<trng::uint128, std::string>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{trng::uint128{0x0, 0x0}, std::string {"0"}},
      data_type{trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}, std::string{"3777777777777777777777777777777777777777777"}},
      data_type{trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}, std::string{"2054255400352543775042143032311251212212673"}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << std::oct << a_y.a;
      REQUIRE(strstr.str() == a_y.y);
    }
    {
      std::stringstream strstr;
      strstr << std::oct << std::showbase << a_y.a;
      REQUIRE(strstr.str() == "0" + a_y.y);
    }
  }
  SECTION("decimal input") {
    using data_type = unary_op_tuple<std::string, trng::uint128>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{std::string {"0"}, trng::uint128{0x0, 0x0}},
      data_type{std::string{"340282366920938463463374607431768211455"}, trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}},
      data_type{std::string{"177508241696884716594379275141355476411"}, trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << "  " << a_y.a << "  ";
      trng::uint128 u;
      strstr >> std::dec >> u;
      REQUIRE(u == a_y.y);
      REQUIRE(strstr.good());
    }
  }
  SECTION("decimal faulty input") {
    SECTION("non digits") {
      std::stringstream strstr;
      strstr << "  non-digits  ";
      trng::uint128 u;
      strstr >> std::dec >> u;
      REQUIRE(strstr.fail());
    }
    SECTION("overflow") {
      std::stringstream strstr;
      strstr << "  340282366920938463463374607431768211456  "; // 2^128
      trng::uint128 u;
      strstr >> std::dec >> u;
      REQUIRE(strstr.fail());
    }
  }
  SECTION("hex input") {
    using data_type = unary_op_tuple<std::string, trng::uint128>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{std::string {"0"}, trng::uint128{0x0, 0x0}},
      data_type{std::string{"ffffffffffffffffffffffffffffffff"}, trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}},
      data_type{std::string{"858ad803aac7fd1118c34c954a2915bb"}, trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << "  " << a_y.a << "  ";
      trng::uint128 u;
      strstr >> std::hex >> u;
      REQUIRE(u == a_y.y);
      REQUIRE(strstr.good());
    }
    {
      std::stringstream strstr;
      strstr << "  0x" << a_y.a << "  ";
      trng::uint128 u;
      strstr >> std::hex >> u;
      REQUIRE(u == a_y.y);
      REQUIRE(strstr.good());
    }
  }
  SECTION("hex faulty input") {
    SECTION("non digits") {
      std::stringstream strstr;
      strstr << "  non-digits  ";
      trng::uint128 u;
      strstr >> std::hex >> u;
      REQUIRE(strstr.fail());
    }
    SECTION("overflow") {
      std::stringstream strstr;
      strstr << "  100000000000000000000000000000000  "; // 2^128
      trng::uint128 u;
      strstr >> std::hex >> u;
      REQUIRE(strstr.fail());
    }
  }
  SECTION("octal input") {
    using data_type = unary_op_tuple<std::string, trng::uint128>;
    // clang-format off
    auto a_y{GENERATE(
      data_type{std::string {"0"}, trng::uint128{0x0, 0x0}},
      data_type{std::string{"3777777777777777777777777777777777777777777"}, trng::uint128{0xffffffffffffffff, 0xffffffffffffffff}},
      data_type{std::string{"2054255400352543775042143032311251212212673"}, trng::uint128{0x858ad803aac7fd11, 0x18c34c954a2915bb}}
    )};
    // clang-format on
    {
      std::stringstream strstr;
      strstr << "  " << a_y.a << "  ";
      trng::uint128 u;
      strstr >> std::oct >> u;
      REQUIRE(u == a_y.y);
      REQUIRE(strstr.good());
    }
  }
  SECTION("octal faulty input") {
    SECTION("non digits") {
      std::stringstream strstr;
      strstr << "  non-digits  ";
      trng::uint128 u;
      strstr >> std::oct >> u;
      REQUIRE(strstr.fail());
    }
    SECTION("overflow") {
      std::stringstream strstr;
      strstr << "  4000000000000000000000000000000000000000000  "; // 2^128
      trng::uint128 u;
      strstr >> std::oct >> u;
      REQUIRE(strstr.fail());
    }
  }
}
