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
#if defined _MSC_VER && __cpluplus <= 201703
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
                   trng::xoshiro256plus) {
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
    SECTION("float") { test_ranges_impl<TestType, float>(); }
    SECTION("double") { test_ranges_impl<TestType, double>(); }
    SECTION("long double") { test_ranges_impl<TestType, long double>(); }
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
        THEN("both engines have equal state") { REQUIRE(r1 == r2); }
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
        THEN("both engines have equal state") { REQUIRE(r1 == r2); }
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
