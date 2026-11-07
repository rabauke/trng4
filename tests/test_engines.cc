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

#include <vector>
#include <string>
#include <sstream>
#include <tuple>
#include <cstdint>
#include <type_traits>
#if defined _MSC_VER && __cplusplus <= 201703
#include <ciso646>
#endif

#include <catch2/catch_all.hpp>

#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/lcg64_count_shift.hpp>
#include <trng/count128_lcg_shift.hpp>
#include <trng/mrg2.hpp>
#include <trng/mrg3.hpp>
#include <trng/mrg3s.hpp>
#include <trng/mrg4.hpp>
#include <trng/mrg5.hpp>
#include <trng/mrg5s.hpp>
#include <trng/yarn2.hpp>
#include <trng/yarn3.hpp>
#include <trng/yarn3s.hpp>
#include <trng/yarn4.hpp>
#include <trng/yarn5.hpp>
#include <trng/yarn5s.hpp>
#include <trng/lagfib2xor.hpp>
#include <trng/lagfib4xor.hpp>
#include <trng/lagfib2plus.hpp>
#include <trng/lagfib4plus.hpp>
#include <trng/mt19937.hpp>
#include <trng/mt19937_64.hpp>
#include <trng/xoshiro256plus.hpp>
#include <trng/minstd.hpp>
#include <trng/mwc64x32.hpp>
#include <trng/mwc128x32.hpp>
#include <trng/mwc128x.hpp>
#include <trng/mwc192x.hpp>


template<typename R>
void advance_engine(R &r1, const long N) {
  for (long i{0}; i < N; ++i)
    r1();
}

template<typename R>
void advance_engine(R &r1, R &r2, const long N) {
  for (long i{0}; i < N; ++i) {
    r1();
    r2();
  }
}


template<typename R>
std::tuple<std::vector<typename R::result_type>, std::vector<typename R::result_type>>
generate_list(R &r1, R &r2, const long N, const long skip1, const long skip2) {
  std::vector<typename R::result_type> v1, v2;
  if (N > 0) {
    v1.reserve(N);
    v2.reserve(N);
  }
  for (long i{0}; i < N; ++i) {
    v1.push_back(r1());
    v2.push_back(r2());
    for (long j{0}; j < skip1; ++j)
      r1();
    for (long j{0}; j < skip2; ++j)
      r2();
  }
  return {v1, v2};
}

template<typename R>
std::tuple<std::vector<typename R::result_type>, std::vector<typename R::result_type>>
generate_list(R &r1, R &r2, const long N) {
  return generate_list(r1, r2, N, 0, 0);
}


template<typename R>
struct generator_min_max {
  typedef typename R::result_type result_type;
  static constexpr result_type min() { return R::min(); }
  static constexpr result_type max() { return R::max(); }
};

template<typename R>
struct generator_min : public generator_min_max<R> {
  using base = generator_min_max<R>;
  typename base::result_type operator()() const { return base::min(); }
};

template<typename R>
struct generator_max : public generator_min_max<R> {
  using base = generator_min_max<R>;
  typename base::result_type operator()() const { return base::max(); }
};


template<typename R, typename T>
void test_ranges_impl() {
  generator_min<R> r_min;
  generator_max<R> r_max;
  using u01xx_min = trng::utility::u01xx_traits<T, 1, generator_min<R>>;
  using u01xx_max = trng::utility::u01xx_traits<T, 1, generator_max<R>>;
  T x_min{}, x_max{};
  x_min = u01xx_min::cc(r_min);
  x_max = u01xx_max::cc(r_max);
  REQUIRE((x_min >= 0 and x_max <= 1));
  x_min = u01xx_min::co(r_min);
  x_max = u01xx_max::co(r_max);
  REQUIRE((x_min >= 0 and x_max < 1));
  x_min = u01xx_min::oc(r_min);
  x_max = u01xx_max::oc(r_max);
  REQUIRE((x_min > 0 and x_max <= 1));
  x_min = u01xx_min::oo(r_min);
  x_max = u01xx_max::oo(r_max);
  REQUIRE((x_min > 0 and x_max < 1));
}


TEMPLATE_TEST_CASE("engines", "",                                            //
                   trng::lcg64, trng::lcg64_shift, trng::lcg64_count_shift,  //
                   trng::count128_lcg_shift,                                 //
                   trng::mrg2, trng::mrg3, trng::mrg3s, trng::mrg4, trng::mrg5,
                   trng::mrg5s,  //
                   trng::yarn2, trng::yarn3, trng::yarn3s, trng::yarn4, trng::yarn5,
                   trng::yarn5s,                                        //
                   trng::lagfib2xor_521_64, trng::lagfib4xor_521_32,    //
                   trng::lagfib2plus_521_32, trng::lagfib4plus_521_64,  //
                   trng::mt19937, trng::mt19937_64,                     //
                   trng::xoshiro256plus,                                //
                   trng::mwc64x32, trng::mwc128x32, trng::mwc128x, trng::mwc192x) {
  SECTION("advance") {
    // two engines with equal state
    TestType r1, r2;
    advance_engine(r1, 271828l);  // advance engine r1
    REQUIRE(r1 != r2);
  }

  SECTION("restore") {
    // three engines with equal state
    TestType r1, r2, r3;
    advance_engine(r1, r2, 271828l);  // advance engines ra and rb
    r3 = r1;                          // backup r1 in r3
    // advance ra such that ra and rb have different states
    advance_engine(r1, 314159l);
    r1 = r3;  // restore ra
    REQUIRE(r1 == r2);
    // when state has been restored correctly, r1 and r2 must yield the same values
    const auto v{generate_list(r1, r2, 32)};
    REQUIRE(std::get<0>(v) == std::get<1>(v));
  }

  SECTION("state_io") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      WHEN("advance 1st engine and restore state from 2nd one") {
        advance_engine(r1, 271828l);  // advance engine r1
        std::stringstream str;
        str << r1;
        str >> r2;
        THEN("engines have same state and generate same values") {
          REQUIRE(r1 == r2);
          // when state has been restored correctly, r1 and r2 must yield the same values
          const auto v{generate_list(r1, r2, 32)};
          REQUIRE(std::get<0>(v) == std::get<1>(v));
        }
      }
      AND_WHEN("restore from empty io-stream") {
        // an empty stream cannot yield a valid engine
        std::stringstream str;
        str >> r1;
        THEN("stream not good and stata unchanged") {
          REQUIRE(not str.good());
          REQUIRE(r1 == r2);
        }
      }
      AND_WHEN("restore from incomplete io-stream") {
        // a stream holding a trimmed status string cannot yield a valid engine
        std::ostringstream strout;
        strout << r1;
        const std::string s{strout.str().substr(0, strout.str().length() / 2)};
        std::stringstream str;
        str << s;
        str >> r1;
        THEN("stream not good and stata unchanged") {
          REQUIRE(not str.good());
          REQUIRE(r1 == r2);
        }
      }
    }
  }

  SECTION("state_wio") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      WHEN("advance 1st engine and restore state from 2nd one") {
        advance_engine(r1, 271828l);  // advance engine r1
        std::wstringstream str;
        str << r1;
        str >> r2;
        THEN("engines have same state and generate same values") {
          REQUIRE(r1 == r2);
          // when state has been restored correctly, r1 and r2 must yield the same values
          const auto v{generate_list(r1, r2, 32)};
          REQUIRE(std::get<0>(v) == std::get<1>(v));
        }
      }
      AND_WHEN("restore from empty io-stream") {
        // an empty stream cannot yield a valid engine
        std::wstringstream str;
        str >> r1;
        THEN("stream not good and stata unchanged") {
          REQUIRE(not str.good());
          REQUIRE(r1 == r2);
        }
      }
      AND_WHEN("restore from incomplete io-stream") {
        // a stream holding a trimmed status string cannot yield a valid engine
        std::wostringstream strout;
        strout << r1;
        const std::wstring s{strout.str().substr(0, strout.str().length() / 2)};
        std::wstringstream str;
        str << s;
        str >> r1;
        THEN("stream not good and stata unchanged") {
          REQUIRE(not str.good());
          REQUIRE(r1 == r2);
        }
      }
    }
  }

  SECTION("ranges") {
    SECTION("float") {
      test_ranges_impl<TestType, float>();
    }
    SECTION("double") {
      test_ranges_impl<TestType, double>();
    }
    SECTION("long double") {
      test_ranges_impl<TestType, long double>();
    }
  }

  SECTION("discard") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      const unsigned long long throw_away =
          GENERATE(2ull * 3 * 5 * 7 * 11 * 13 * 17 * 19 + 0x10000000ull,
                   2ull * 3 * 5 * 7 * 11 * 13 * 17 * 19);
      WHEN("discard one, advance other") {
        r1.discard(throw_away);
        for (unsigned long long i{0}; i < throw_away; ++i)
          r2();
        THEN("both engines have equal state") {
          REQUIRE(r1 == r2);
        }
      }
    }
  }
}


TEMPLATE_TEST_CASE("parallel engines", "",                                   //
                   trng::lcg64, trng::lcg64_shift, trng::lcg64_count_shift,  //
                   trng::count128_lcg_shift,                                 //
                   trng::mrg2, trng::mrg3, trng::mrg3s, trng::mrg4, trng::mrg5,
                   trng::mrg5s,  //
                   trng::yarn2, trng::yarn3, trng::yarn3s, trng::yarn4, trng::yarn5,
                   trng::yarn5s) {
  SECTION("jump2") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      const long i{GENERATE(range(0l, 20l))};
      const long n{1l << i};
      WHEN("jump ahead one, advance other") {
        r1.jump2(i);
        for (long j{0l}; j < n; ++j)
          r2();
        THEN("both engines have equal state") {
          REQUIRE(r1 == r2);
        }
      }
    }
  }

  SECTION("split") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      const long i{GENERATE(range(2l, 21l))};
      const long j{GENERATE_COPY(range(0l, i))};
      WHEN("split one, advance other") {
        for (long k{0l}; k < j; ++k)
          r1();
        r2.split(i, j);  // split into i streams, pick stream j
        THEN(
            "split engine yields the same values as when values are skipped in the other "
            "engine") {
          const auto v{generate_list(r1, r2, 32, i - 1, 0)};
          REQUIRE(std::get<0>(v) == std::get<1>(v));
        }
      }
    }
  }
}


TEMPLATE_TEST_CASE("block splitting engines", "", trng::xoshiro256plus, trng::mwc64x32,
                   trng::mwc128x32, trng::mwc128x, trng::mwc192x) {
  SECTION("jump2") {
    // two engines with equal state
    GIVEN("two engines with equal state") {
      TestType r1, r2;
      const long i{GENERATE(range(0l, 20l))};
      const long n{1l << i};
      WHEN("jump ahead one, advance other") {
        r1.jump2(i);
        for (long j{0l}; j < n; ++j)
          r2();
        THEN("both engines have equal state") {
          REQUIRE(r1 == r2);
        }
      }
    }
  }
}


// Reference values have been generated with the C code from https://github.com/vigna/MWCPRNG.
// The raw state is loaded via the stream interface to test reference states without seeding.
template<typename R>
R engine_from_state(const std::string &state) {
  R r;
  std::stringstream str;
  str << '[' << R::name() << ' ' << state << ']';
  str >> r;
  REQUIRE(static_cast<bool>(str));
  return r;
}


template<typename R>
std::string engine_state(const R &r) {
  std::stringstream str;
  str << r;
  std::string s{str.str()};
  // strip leading "[name " and trailing "]"
  s = s.substr(std::string(R::name()).size() + 2, std::string::npos);
  s.pop_back();
  return s;
}


TEST_CASE("mwc64x32 reference vectors") {
  using R = trng::mwc64x32;
  const std::string state0{"(2309737967 19088743)"};  // x=0x89abcdef, c=0x01234567
  SECTION("output") {
    R r{engine_from_state<R>(state0)};
    const R::result_type expected[]{0x88888888u, 0x14b0cd91u, 0x8c2e9bb8u, 0x1d1a2201u,
                                    0x30ccae9eu};
    for (const R::result_type e : expected)
      REQUIRE(r() == e);
    REQUIRE(engine_state(r) == "(3284445065 1022492261)");
  }
  SECTION("jump2") {
    struct {
      unsigned e;
      const char *state;
      R::result_type next;
    } const ref[]{{0, "(2635280367 2309218942)", 0x14b0cd91u},
                  {1, "(287620358 2634688190)", 0x8c2e9bb8u},
                  {10, "(884098881 3050340339)", 0x816238b2u},
                  {32, "(3228847466 2961994600)", 0x70f83e02u},
                  {48, "(2005114421 4045125995)", 0x8698275eu},
                  {63, "(3657244213 1515460187)", 0x83a9086eu}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump2(t.e);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
  SECTION("jump") {
    struct {
      unsigned long long n;
      const char *state;
      R::result_type next;
    } const ref[]{{63, "(3032529233 927382902)", 0x83860427u},
                  {64, "(1747851502 3031847789)", 0xdc984b83u},
                  {65, "(2552537085 1747458740)", 0xf00c8549u},
                  {1000, "(775641254 3517362959)", 0xff9deba9u},
                  {123456789012345ull, "(2123935596 177324147)", 0x74096f1fu}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump(t.n);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
}


TEST_CASE("mwc128x32 reference vectors") {
  using R = trng::mwc128x32;
  const std::string state0{"(2309737967 324508639 4275878552 19088743)"};
  SECTION("output") {
    R r{engine_from_state<R>(state0)};
    const R::result_type expected[]{0xffffffffu, 0x66b58255u, 0x7631ae4bu, 0x63bef72au,
                                    0xe7bec082u};
    for (const R::result_type e : expected)
      REQUIRE(r() == e);
    REQUIRE(engine_state(r) == "(2635544085 155922789 1979317755 1699597525)");
  }
  SECTION("jump2") {
    struct {
      unsigned e;
      const char *state;
      R::result_type next;
    } const ref[]{{0u, "(324508639 4275878552 4012230305 2307933428)", 0x66b58255u},
                  {1, "(4275878552 4012230305 1700926414 324255109)", 0x7631ae4bu},
                  {10, "(3954778911 974194992 4033717946 3478049772)", 0x3f237556u},
                  {32, "(3442286283 3536530296 1474472042 464651567)", 0x4c50ad45u},
                  {64, "(244646757 1777918563 739090470 3503274784)", 0xfcc21b06u},
                  {96, "(3292317958 4007266617 943325353 722075108)", 0x1333fd4du},
                  {100, "(506206589 782810952 1654225643 4088951221)", 0x91210b5eu}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump2(t.e);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
  SECTION("jump") {
    struct {
      unsigned long long n;
      const char *state;
      R::result_type next;
    } const ref[]{
        {63, "(847885179 3715962022 2476216921 3652926754)", 0x4a232f7bu},
        {64, "(3715962022 2476216921 1117963300 847222749)", 0x70dd53f9u},
        {65, "(2476216921 1117963300 3977386497 3713058837)", 0x3042d214u},
        {1000, "(424111210 3721378245 3644612862 2103390486)", 0xa46351e8u},
        {123456789012345ull, "(4075275988 3027750538 3476453331 4211364144)", 0x343226e3u}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump(t.n);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
}


TEST_CASE("mwc128x reference vectors") {
  using R = trng::mwc128x;
  const std::string state0{"(81985529216486895 1393753991812647630)"};
  SECTION("output") {
    R r{engine_from_state<R>(state0)};
    const R::result_type expected[]{0x1274deb88bed4721u, 0xcfae465e857bb779u,
                                    0xec6f86cb0882e858u, 0x5460b4c476d45ae8u,
                                    0x47f15b994cb2364eu};
    for (const R::result_type e : expected)
      REQUIRE(r() == e);
    REQUIRE(engine_state(r) == "(13671451619669037541 3608578230611311617)");
  }
  SECTION("jump2") {
    struct {
      unsigned e;
      const char *state;
      R::result_type next;
    } const ref[]{{0u, "(14883508920264627904 81802424913797561)", 0xcfae465e857bb779u},
                  {1, "(2484093949676386361 14850268486880104545)", 0xec6f86cb0882e858u},
                  {10, "(1501402912552787797 13596553788164794155)", 0xa866a45f7a7ff87eu},
                  {64, "(17220119740164636705 14890220014685849057)", 0x205e979cccf785c0u},
                  {96, "(16050319296952212815 16125258893842661718)", 0x01765ddb94dc5c19u},
                  {100, "(8860938554697113021 3614744106532962313)", 0x48d24637e69d5db4u}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump2(t.e);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
  SECTION("jump") {
    struct {
      unsigned long long n;
      const char *state;
      R::result_type next;
    } const ref[]{{63, "(3370160081764681560 8895604331013298097)", 0x55b6bfbf19b050e9u},
                  {64, "(14269024193373299073 3362633255779409339)", 0xe8afb8c611ddfc3au},
                  {65, "(16730373316449978345 14237156133851650724)", 0x2dbab3f1a4ac9d4du},
                  {1000, "(13012910374931284489 12870968565340708924)", 0x0609f8bde39c6235u},
                  {123456789012345ull, "(10952109836232077281 4212875018435523883)",
                   0xad8a965ca55bdacau}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump(t.n);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
}


TEST_CASE("mwc192x reference vectors") {
  using R = trng::mwc192x;
  const std::string state0{"(81985529216486895 18364758544493064720 1393753991812647630)"};
  SECTION("output") {
    R r{engine_from_state<R>(state0)};
    const R::result_type expected[]{0xed8b21477412b8deu, 0xc025294bbfbdd881u,
                                    0xaa0f361e7b30715du, 0x3652fcc212c0ab04u,
                                    0xc0a8dde1a134ab38u};
    for (const R::result_type e : expected)
      REQUIRE(r() == e);
    REQUIRE(engine_state(r) ==
            "(10611315881025267658 7360190193432640826 17737513442500101382)");
  }
  SECTION("jump2") {
    struct {
      unsigned e;
      const char *state;
      R::result_type next;
    } const ref[]{{0, "(18364758544493064720 13909176411246296427 81763031380272618)",
                   0xc025294bbfbdd881u},
                  {1, "(13909176411246296427 6063146698804387610 18314919029181066823)",
                   0xaa0f361e7b30715du},
                  {10, "(9863171486552869619 10278578489538821469 7287101682763666229)",
                   0xeb8424f2bba92e68u},
                  {64, "(18112313878986340860 1271982935409551474 12002627988787756372)",
                   0xb7371b2dbbab1526u},
                  {96, "(18191567408486071481 15233050496709431133 17234047554916267038)",
                   0x3c4d08470d0d3b43u},
                  {144, "(16474492249774699218 13783717589596303460 7484804607378852770)",
                   0xd896c43c0836fbc6u},
                  {150, "(8456223794033853695 2192494831446734257 14119284508329342731)",
                   0xdd9c881ad1e642bau}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump2(t.e);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
  SECTION("jump") {
    struct {
      unsigned long long n;
      const char *state;
      R::result_type next;
    } const ref[]{
        {63, "(8176644946430835372 15873698107482791753 10538621308507988564)",
         0x4e0a006df3b67d1du},
        {64, "(15873698107482791753 1136558464365281944 8154454618144128301)",
         0x7eefaca3ad72ffb5u},
        {65, "(1136558464365281944 1363487952914295992 15830618999311061693)",
         0xc95dbd0abc877e05u},
        {1000, "(14637051160347023684 13691764478331869955 1194771707803656226)",
         0xae9645c9801abb21u},
        {123456789012345ull, "(3717131272906053156 10000394222409139780 8402368892279984879)",
         0xfe53bc72f0b4d4abu}};
    for (const auto &t : ref) {
      R r{engine_from_state<R>(state0)};
      r.jump(t.n);
      REQUIRE(engine_state(r) == t.state);
      REQUIRE(r() == t.next);
    }
  }
}


template<typename R>
struct mwc_traits;

template<>
struct mwc_traits<trng::mwc64x32> {
  static constexpr unsigned long long A{0xfff145f8ull};
  static void seed_raw(trng::mwc64x32 &r, std::uint64_t c) {
    r.seed(0x89abcdefu, static_cast<std::uint32_t>(c));
  }
  static const char *noncanonical_state() { return "(4294967295 4294967295)"; }
};

template<>
struct mwc_traits<trng::mwc128x32> {
  static constexpr unsigned long long A{0xffcccc66ull};
  static void seed_raw(trng::mwc128x32 &r, std::uint64_t c) {
    r.seed(0x89abcdefu, 0x13579bdfu, 0xfedcba98u, static_cast<std::uint32_t>(c));
  }
  static const char *noncanonical_state() {
    return "(4294967295 4294967295 4294967295 4294967295)";
  }
};

template<>
struct mwc_traits<trng::mwc128x> {
  static constexpr unsigned long long A{0xff6da236274e252eull};
  static void seed_raw(trng::mwc128x &r, std::uint64_t c) { r.seed(0x0123456789abcdefull, c); }
  static const char *noncanonical_state() {
    return "(18446744073709551615 18446744073709551615)";
  }
};

template<>
struct mwc_traits<trng::mwc192x> {
  static constexpr unsigned long long A{0xff4e24dd418f5433ull};
  static void seed_raw(trng::mwc192x &r, std::uint64_t c) {
    r.seed(0x0123456789abcdefull, 0xfedcba9876543210ull, c);
  }
  static const char *noncanonical_state() {
    return "(18446744073709551615 18446744073709551615 18446744073709551615)";
  }
};


// the carry is the last number in the textual representation of the engine state
template<typename R>
unsigned long long engine_carry(const R &r) {
  const std::string s{engine_state(r)};
  const std::size_t pos{s.find_last_of(' ')};
  return std::stoull(s.substr(pos + 1));
}

template<typename R>
bool carry_valid(const R &r) {
  const unsigned long long c{engine_carry(r)};
  return c > 0 and c < mwc_traits<R>::A - 1;
}


TEMPLATE_TEST_CASE("mwc engines carry and seeding", "", trng::mwc64x32, trng::mwc128x32,
                   trng::mwc128x, trng::mwc192x) {
  using R = TestType;
  const unsigned long long A{mwc_traits<R>::A};
  const unsigned long long cmax{static_cast<unsigned long long>(R::max())};

  SECTION("default state is valid") {
    REQUIRE(carry_valid(R()));
  }

  SECTION("raw seed preserves valid carries and maps invalid values") {
    const unsigned long long cs[]{0, 1, 2, A - 3, A - 2, A - 1, A, A + 1, cmax - 1, cmax};
    for (const unsigned long long c : cs) {
      R r;
      mwc_traits<R>::seed_raw(r, c);
      REQUIRE(carry_valid(r));
      const unsigned long long expected{c > 0 && c < A - 1 ? c : c % (A - 2) + 1};
      REQUIRE(engine_carry(r) == expected);
      // raw seeding is deterministic
      R r2;
      mwc_traits<R>::seed_raw(r2, c);
      REQUIRE(r == r2);
    }
  }

  SECTION("zero carry after stepping is restorable and jumpable") {
    R r;
    r();
    REQUIRE(engine_carry(r) == 0);

    std::stringstream state;
    state << r;
    R restored;
    state >> restored;
    REQUIRE(state);
    REQUIRE(restored == r);

    R jumped = restored;
    R stepped = restored;
    jumped.jump(64);
    for (int i{0}; i < 64; ++i)
      stepped();
    REQUIRE(jumped == stepped);
  }

  SECTION("noncanonical stream states are rejected without assignment") {
    R r;
    R original = r;
    std::stringstream state;
    state << '[' << R::name() << ' ' << mwc_traits<R>::noncanonical_state() << ']';
    state >> r;
    REQUIRE(state.fail());
    REQUIRE(r == original);
  }

  SECTION("seeded states are valid") {
    for (unsigned long s{0}; s < 2000; ++s) {
      R r(s);
      REQUIRE(carry_valid(r));
    }
    trng::minstd g(12345);
    for (int i{0}; i < 500; ++i) {
      R r(g);
      REQUIRE(carry_valid(r));
    }
    R r(5ul);
    r.seed();
    REQUIRE(r == R());
  }

  SECTION("different seeds give different states") {
    REQUIRE(R(1ul) != R(2ul));
    REQUIRE(R(1ul) == R(1ul));
  }
}


TEMPLATE_TEST_CASE("mwc engines jump", "", trng::mwc64x32, trng::mwc128x32, trng::mwc128x,
                   trng::mwc192x) {
  using R = TestType;
  SECTION("jump equals stepping around the threshold") {
    const unsigned long long ns[]{0, 1, 2, 31, 62, 63, 64, 65, 66, 100, 127, 128, 1000, 4097};
    for (const unsigned long long n : ns) {
      R r1(42ul), r2(42ul);
      r1.jump(n);
      for (unsigned long long i{0}; i < n; ++i)
        r2();
      REQUIRE(r1 == r2);
      REQUIRE(r1() == r2());
    }
  }

  SECTION("jump of sum equals consecutive jumps") {
    R r1(7ul), r2(7ul);
    r1.jump(1000000007ull);
    r2.jump(1000000000ull);
    r2.jump(7);
    REQUIRE(r1 == r2);
    r1.jump(123456789012345ull);
    r2.jump(100000000000000ull);
    r2.jump(23456789012345ull);
    REQUIRE(r1 == r2);
  }

  SECTION("jump2 equals jump of power of two") {
    for (unsigned e{0}; e < 40; e += 3) {
      R r1(9ul), r2(9ul);
      r1.jump2(e);
      r2.jump(1ull << e);
      REQUIRE(r1 == r2);
    }
  }

  SECTION("jumping keeps the carry valid") {
    R r(3ul);
    for (unsigned e{0}; e < 200; e += 7) {
      r.jump2(e);
      REQUIRE(carry_valid(r));
    }
  }
}
