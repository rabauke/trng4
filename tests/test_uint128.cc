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

#undef __SIZEOF_INT128__
#include <trng/uint128.hpp>

struct tuple {
  trng::uint128 a;
  trng::uint128 b;
  trng::uint128 y;
};


TEST_CASE("uint128 math") {
  SECTION("arithmetics") {
    SECTION("plus") {
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0xd091bb5c22ae9ef6, 0xe7e1faeed5c31f79}, trng::uint128{0xa7de9f4ccc450cba, 0x924668f5c7dc380}, trng::uint128{0x78705aa8eef3abb0, 0xf106617e3240e2f9}},
          tuple{trng::uint128{0x2082352cf807b7df, 0xe9d300053895afe1}, trng::uint128{0xd96089c53640ac4c, 0xef1a2e6dae6d9426}, trng::uint128{0xf9e2bef22e48642c, 0xd8ed2e72e7034407}},
          tuple{trng::uint128{0xa1e24bba4ee4092b, 0x18f868638c16a625}, trng::uint128{0xadc1965b6613ba46, 0xc1fb41c2bd9b0ecd}, trng::uint128{0x4fa3e215b4f7c371, 0xdaf3aa2649b1b4f2}},
          tuple{trng::uint128{0x474ba8c43039cd1a, 0x8c006d5ffe2d7810}, trng::uint128{0xbe3dedfc7989c8ee, 0x6468fd6e6c0df032}, trng::uint128{0x58996c0a9c39608, 0xf0696ace6a3b6842}},
          tuple{trng::uint128{0xf51f2ae7ff1816e4, 0xf702ef59f7badafa}, trng::uint128{0xa7cd66342c826d8b, 0x2bd2e4124d4a2dbe}, trng::uint128{0x9cec911c2b9a8470, 0x22d5d36c450508b8}},
          tuple{trng::uint128{0x285954a1b9d09511, 0xf878c4b3fb2a0137}, trng::uint128{0xb4bf6fa7cc1a8959, 0x826328251097330}, trng::uint128{0xdd18c44985eb1e6b, 0x9ef7364c337467}},
          tuple{trng::uint128{0xf508e4aa1c1fe652, 0x7c419418cc50aa59}, trng::uint128{0x46e46cb0df577ec2, 0xbd1e364262c5564}, trng::uint128{0x3bed515afb776514, 0x8813777cf27cffbd}},
          tuple{trng::uint128{0xccdf2e5c4c0a1f3b, 0x2452a9dc01397d8d}, trng::uint128{0x18dda0c9fe7b45d9, 0xd2ce21c9d268409a}, trng::uint128{0xe5bccf264a856514, 0xf720cba5d3a1be27}},
          tuple{trng::uint128{0x6bf88c311cca797a, 0xea6da4aea3c78807}, trng::uint128{0xb1e049e1200bfa47, 0x512d6e73c3851eee}, trng::uint128{0x1dd8d6123cd673c2, 0x3b9b1322674ca6f5}},
          tuple{trng::uint128{0xcace1969e0e0d4ad, 0xf5a14bab80f00988}, trng::uint128{0xf341c0817d973e48, 0x8d17554a9e20d28}, trng::uint128{0xbe0fd9eb5e7812f5, 0xfe72c1002ad216b0}}
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
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0x70518ce6203ac303, 0x61add0ab35d0430c}, trng::uint128{0xc05309bed23d2d63, 0x414de9c5d2229f23}, trng::uint128{0xaffe83274dfd95a0, 0x205fe6e563ada3e9}},
          tuple{trng::uint128{0xc3f8e8920d1c8509, 0xcb92388e095436bf}, trng::uint128{0x818666a3f0a8b109, 0xb2f6b12769a48341}, trng::uint128{0x427281ee1c73d400, 0x189b87669fafb37e}},
          tuple{trng::uint128{0x2fd6e20868a29af9, 0x7d61330b753ec6fc}, trng::uint128{0xe4123c566c548c8f, 0xf5941f6194b993aa}, trng::uint128{0x4bc4a5b1fc4e0e69, 0x87cd13a9e0853352}},
          tuple{trng::uint128{0x7211efea7cd15133, 0xa574c4ffcb41f198}, trng::uint128{0x8c1651342876763c, 0x237ce42ec300d11b}, trng::uint128{0xe5fb9eb6545adaf7, 0x81f7e0d10841207d}},
          tuple{trng::uint128{0xb598eef6ebbe7347, 0xc1332568ceba5a70}, trng::uint128{0x263821ca3aeb8202, 0x41ec0f84cf4ac36d}, trng::uint128{0x8f60cd2cb0d2f145, 0x7f4715e3ff6f9703}},
          tuple{trng::uint128{0x46a99459b4ad9f11, 0xae00feaa00b8b573}, trng::uint128{0xd7393ee6fd0fc06a, 0x4118a30a551b54a4}, trng::uint128{0x6f705572b79ddea7, 0x6ce85b9fab9d60cf}},
          tuple{trng::uint128{0xa7b480b6b5f0b06c, 0x29a0ec27a4daa010}, trng::uint128{0xd074f86f4cc1c54a, 0x3e57a70303774cda}, trng::uint128{0xd73f8847692eeb21, 0xeb494524a1635336}},
          tuple{trng::uint128{0x1e76a1c574be9133, 0x7f94c950c61f6ed6}, trng::uint128{0xede43895379ce627, 0x59988939e8490ddc}, trng::uint128{0x309269303d21ab0c, 0x25fc4016ddd660fa}},
          tuple{trng::uint128{0xf5b1c7a192e195f8, 0x572384d4e0732c88}, trng::uint128{0x325410e1d9352f6a, 0x4047080af47c081d}, trng::uint128{0xc35db6bfb9ac668e, 0x16dc7cc9ebf7246b}},
          tuple{trng::uint128{0x95d41b68cee496c3, 0x394bbd52048cd47c}, trng::uint128{0x9db51a85c765d71f, 0x79297527fcca2773}, trng::uint128{0xf81f00e3077ebfa3, 0xc022482a07c2ad09}}
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
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0x70518ce6203ac303, 0x61add0ab35d0430c}, trng::uint128{0xc05309bed23d2d63, 0x414de9c5d2229f23}, trng::uint128{0xaffe83274dfd95a0, 0x205fe6e563ada3e9}},
          tuple{trng::uint128{0xc3f8e8920d1c8509, 0xcb92388e095436bf}, trng::uint128{0x818666a3f0a8b109, 0xb2f6b12769a48341}, trng::uint128{0x427281ee1c73d400, 0x189b87669fafb37e}},
          tuple{trng::uint128{0x2fd6e20868a29af9, 0x7d61330b753ec6fc}, trng::uint128{0xe4123c566c548c8f, 0xf5941f6194b993aa}, trng::uint128{0x4bc4a5b1fc4e0e69, 0x87cd13a9e0853352}},
          tuple{trng::uint128{0x7211efea7cd15133, 0xa574c4ffcb41f198}, trng::uint128{0x8c1651342876763c, 0x237ce42ec300d11b}, trng::uint128{0xe5fb9eb6545adaf7, 0x81f7e0d10841207d}},
          tuple{trng::uint128{0xb598eef6ebbe7347, 0xc1332568ceba5a70}, trng::uint128{0x263821ca3aeb8202, 0x41ec0f84cf4ac36d}, trng::uint128{0x8f60cd2cb0d2f145, 0x7f4715e3ff6f9703}},
          tuple{trng::uint128{0x46a99459b4ad9f11, 0xae00feaa00b8b573}, trng::uint128{0xd7393ee6fd0fc06a, 0x4118a30a551b54a4}, trng::uint128{0x6f705572b79ddea7, 0x6ce85b9fab9d60cf}},
          tuple{trng::uint128{0xa7b480b6b5f0b06c, 0x29a0ec27a4daa010}, trng::uint128{0xd074f86f4cc1c54a, 0x3e57a70303774cda}, trng::uint128{0xd73f8847692eeb21, 0xeb494524a1635336}},
          tuple{trng::uint128{0x1e76a1c574be9133, 0x7f94c950c61f6ed6}, trng::uint128{0xede43895379ce627, 0x59988939e8490ddc}, trng::uint128{0x309269303d21ab0c, 0x25fc4016ddd660fa}},
          tuple{trng::uint128{0xf5b1c7a192e195f8, 0x572384d4e0732c88}, trng::uint128{0x325410e1d9352f6a, 0x4047080af47c081d}, trng::uint128{0xc35db6bfb9ac668e, 0x16dc7cc9ebf7246b}},
          tuple{trng::uint128{0x95d41b68cee496c3, 0x394bbd52048cd47c}, trng::uint128{0x9db51a85c765d71f, 0x79297527fcca2773}, trng::uint128{0xf81f00e3077ebfa3, 0xc022482a07c2ad09}}
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
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0x5a065b97114dee4f, 0xd4b12f5fcb29360a}, trng::uint128{0x2984c787ed702bbe, 0xcb563b4d6fa56696}, trng::uint128{0x24763a4b1d26603, 0xfabc4c13a01fa5dc}},
          tuple{trng::uint128{0x95d3de16983162a8, 0x8cbaafb3bb98b27f}, trng::uint128{0x4fabc9eddcd87a48, 0x874df2959ecfe9f0}, trng::uint128{0xe2536339b924c138, 0xadb2ecb904dee10}},
          tuple{trng::uint128{0xeacd3439b1fac842, 0x492cbef1ae08ab78}, trng::uint128{0x2a67f49f1e9aa4e1, 0x9a1b7d0878d22934}, trng::uint128{0xa6f5e3f20c286b74, 0xe5274533a5a90c60}},
          tuple{trng::uint128{0xc1d7dfd0646f1d40, 0xc0f463c48fc23a81}, trng::uint128{0x435216025718a361, 0xa771ba4487a3b97c}, trng::uint128{0xc73f32843d47a76, 0x887a98e2457e8f7c}},
          tuple{trng::uint128{0x6164e6233543f2bc, 0x915cc2538701d0df}, trng::uint128{0xb0705c82b7526048, 0xbf86dcd7fd066ea4}, trng::uint128{0x5696e5022a86fc29, 0xbd0e66458d23a0dc}},
          tuple{trng::uint128{0x136b2fdd677a359e, 0xdcfacd05a4ea31e}, trng::uint128{0x7356b1bbb872426d, 0x1575515de99eadb3}, trng::uint128{0x89d1f2029c7474aa, 0xfdabb19b43bb53fa}},
          tuple{trng::uint128{0x87e2593597c34e42, 0xc77780f05b396fba}, trng::uint128{0x3a9e3c0f8168599c, 0xe9d07a328eeab382}, trng::uint128{0x1f580c875017eb41, 0x49bb3ea4c84dca74}},
          tuple{trng::uint128{0xef1b52e6f7080941, 0x2141888b278946b0}, trng::uint128{0x27023ee880d10fac, 0xd368bdc27664b5a7}, trng::uint128{0x659fb28888348175, 0x9896af4f96478cd0}},
          tuple{trng::uint128{0x919e6d646518b459, 0x7829fc226325d30e}, trng::uint128{0x89d0cf468bed7368, 0xff02af497294e430}, trng::uint128{0x947f48a83fc8d24, 0xb9a1d19887280aa0}},
          tuple{trng::uint128{0x30c0399ba19b463, 0x564dab7563794f97}, trng::uint128{0x14034fbbdabd4cc4, 0x71535cf89aaeea20}, trng::uint128{0xa6f23c77915994c7, 0x5b23b036408bf8e0}}
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
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0x1b4d989d7fa09780, 0xf63ef3d2fadc6788}, trng::uint128{0xda603f4888a, 0xfd7149f3f014d704}, trng::uint128{0x0, 0x2001d}},
          tuple{trng::uint128{0x12fb56808c904fa, 0xc660883ffa1cce2a}, trng::uint128{0x59d803fb331, 0x5c9dc83645273235}, trng::uint128{0x0, 0x3616}},
          tuple{trng::uint128{0xd13ac8b85cf9c9b3, 0xde62c6bdadf500ad}, trng::uint128{0x66dce43761e, 0x924413728d90c32}, trng::uint128{0x0, 0x208b86}},
          tuple{trng::uint128{0x159d967e58a2c06c, 0x665827cbdb1aa208}, trng::uint128{0x3d6b25290a5, 0xc50941f91d46440c}, trng::uint128{0x0, 0x5a18a}},
          tuple{trng::uint128{0x4286ddf10b8905b4, 0xccd149a4a8fd9757}, trng::uint128{0x2f1494b9813, 0x8ac023d6d8755d29}, trng::uint128{0x0, 0x169bd7}},
          tuple{trng::uint128{0x6e7122f0bffc21b1, 0xe9203368220c0724}, trng::uint128{0x6ad203b6fb5, 0x234e00cb62703d2c}, trng::uint128{0x0, 0x108ade}},
          tuple{trng::uint128{0x2e8d86cbfb7bfb5e, 0x438896871869325b}, trng::uint128{0xe718672d4d4, 0x48df1f1dd92d70c4}, trng::uint128{0x0, 0x3391d}},
          tuple{trng::uint128{0x25420afc485d46db, 0x22d56381cd572d60}, trng::uint128{0x7da944f13f7, 0x2f0c37d3fa9308b4}, trng::uint128{0x0, 0x4be72}},
          tuple{trng::uint128{0xde89ef2b13dac708, 0x9467851da09c428d}, trng::uint128{0x5674c771e34, 0x20514e669edd5580}, trng::uint128{0x0, 0x292f22}},
          tuple{trng::uint128{0x8cc3a36c0212714e, 0x251bc1f6ae274af0}, trng::uint128{0x5e86f504088, 0x449f41c77c8a029b}, trng::uint128{0x0, 0x17d384}}
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
      // clang-format off
      auto a_b_y{GENERATE(
          tuple{trng::uint128{0xc7be9961e09aebe7, 0x63c5ecb935d657e1}, trng::uint128{0x8c08c64db7e, 0xda5894bdbae3349a}, trng::uint128{0xc2cadb34c7, 0xa7165e947c867b45}},
          tuple{trng::uint128{0x3ddf7ae445d3249c, 0x6766c9407e10ad9c}, trng::uint128{0xbea4300c9f9, 0x97956305f5b0ccfd}, trng::uint128{0xf2d6010901, 0xcf0fdf9b64a8c8b3}},
          tuple{trng::uint128{0x18b13e7f39320ca2, 0x21c900787d84661c}, trng::uint128{0xafd083db500, 0x85e02efa964462b}, trng::uint128{0x8c377809ad4, 0xc895ee46eeab73db}},
          tuple{trng::uint128{0xf12a3a21f4772b41, 0xf4c53baba6e76b3b}, trng::uint128{0x5e5501a5580, 0xbc7bba02889add0d}, trng::uint128{0x3d58453c760, 0x80447095dc3313e7}},
          tuple{trng::uint128{0x9340ded2c1ebc21f, 0xf4db6546f6c42e3}, trng::uint128{0xc7bcfc3cab5, 0x6454f14c4a882612}, trng::uint128{0x916847bdf14, 0x2e4a4b08b9435e4d}},
          tuple{trng::uint128{0x3c1a89438d899f74, 0x5a6899ba0d9b6827}, trng::uint128{0xedec3d2f756, 0xd88c6951b284db}, trng::uint128{0x755ef1382e3, 0x16c9bda2f6fcd7e4}},
          tuple{trng::uint128{0xd239c5cb5290106a, 0x3f17adb67acdc26}, trng::uint128{0x7c9e604de91, 0xac9746f946da27d1}, trng::uint128{0x1c37e22c5a8, 0x2f6f6bbf58a14a95}},
          tuple{trng::uint128{0xb039b90e88f1afa, 0x2b42ee31cf239e4d}, trng::uint128{0x72606b78c17, 0xb72574e6ceb4e1ba}, trng::uint128{0x2e3881468e9, 0xd33f78e29a025c61}},
          tuple{trng::uint128{0xa62c6e93421fae11, 0xbb522891213d9f32}, trng::uint128{0x822d9f9a64c, 0x24df682c1eccbdec}, trng::uint128{0x719c7aa2f56, 0x2ba4523e3b1d00ae}},
          tuple{trng::uint128{0xa5d2adeb7e4aab21, 0x736fbc7560e56773}, trng::uint128{0xd1500cb65c6, 0x15e22cb7a1247eac}, trng::uint128{0x2c04df82740, 0x735402bdb3cc9cd7}}
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
