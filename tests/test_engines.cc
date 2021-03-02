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

#include <vector>
#include <string>
#include <sstream>
#include <tuple>
#include <ciso646>

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 30
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/lcg64_count_shift.hpp>
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

#include "type_names.hpp"

using engines =
    boost::mpl::list<trng::lcg64, trng::lcg64_shift, trng::lcg64_count_shift,  //
                     trng::mrg2, trng::mrg3, trng::mrg3s, trng::mrg4, trng::mrg5,
                     trng::mrg5s,  //
                     trng::yarn2, trng::yarn3, trng::yarn3s, trng::yarn4, trng::yarn5,
                     trng::yarn5s,                                        //
                     trng::lagfib2xor_521_64, trng::lagfib4xor_521_32,    //
                     trng::lagfib2plus_521_32, trng::lagfib4plus_521_64,  //
                     trng::mt19937, trng::mt19937_64,                     //
                     trng::xoshiro256plus>;

using parallel_engines =
    boost::mpl::list<trng::lcg64, trng::lcg64_shift, trng::lcg64_count_shift,  //
                     trng::mrg2, trng::mrg3, trng::mrg3s, trng::mrg4, trng::mrg5,
                     trng::mrg5s,  //
                     trng::yarn2, trng::yarn3, trng::yarn3s, trng::yarn4, trng::yarn5,
                     trng::yarn5s>;

using floats = boost::mpl::list<float, double, long double>;


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

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_engines)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_advance)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_advance, R, engines) {
  // two engines with equal state
  R r1, r2;
  advance_engine(r1, 271828l);  // advance engine r1
  BOOST_TEST(r1 != r2, "engines have different state");
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_restore)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_restore, R, engines) {
  // three engines with equal state
  R r1, r2, r3;
  advance_engine(r1, r2, 271828l);  // advance engines ra and rb
  r3 = r1;                          // backup r1 in r3
  // advance ra such that ra and rb have different states
  advance_engine(r1, 314159l);
  r1 = r3;  // restore ra
  BOOST_TEST(r1 == r2, "engine state restored successfully");
  // when state has been restored correctly, r1 and r2 must yield the same values
  const auto v{generate_list(r1, r2, 32)};
  BOOST_TEST(std::get<0>(v) == std::get<1>(v), "engines yield same values after restore");
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_status_io)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_status_io, R, engines) {
  // two engines with equal state
  R r1, r2;
  {
    advance_engine(r1, 271828l);  // advance engine r1
    std::stringstream str;
    str << r1;
    str >> r2;
    BOOST_TEST(r1 == r2, "engine state restored from stream successfully");
    // when state has been restored correctly, r1 and r2 must yield the same values
    const auto v{generate_list(r1, r2, 32)};
    BOOST_TEST(std::get<0>(v) == std::get<1>(v),
               "engines yield same values after restore from stream");
  }
  {
    // an empty stream cannot yield a valid engine
    std::stringstream str;
    str >> r1;
    BOOST_TEST(not str.good(), "stream is not good");
    BOOST_TEST(r1 == r2, "engine not changed");
  }
  {
    // a stream holding a trimmed status string cannot yield a valid engine
    std::ostringstream strout;
    strout << r1;
    const std::string s{strout.str().substr(0, 4)};
    std::stringstream str;
    str << s;
    str >> r1;
    BOOST_TEST(not str.good(), "stream is not good");
    BOOST_TEST(r1 == r2, "engine not changed");
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_status_wio, R, engines) {
  // two engines with equal state
  R r1, r2;
  {
    advance_engine(r1, 271828l);  // advance engine r1
    std::wstringstream str;
    str << r1;
    str >> r2;
    BOOST_TEST(r1 == r2, "engine state restored from stream successfully");
    // when state has been restored correctly, r1 and r2 must yield the same values
    const auto v{generate_list(r1, r2, 32)};
    BOOST_TEST(std::get<0>(v) == std::get<1>(v),
               "engines yield same values after restore from stream");
  }
  {
    // an empty stream cannot yield a valid engine
    std::wstringstream str;
    str >> r1;
    BOOST_TEST(not str.good(), "stream is not good");
    BOOST_TEST(r1 == r2, "engine not changed");
  }
  {
    // a stream holding a trimmed status string cannot yield a valid engine
    std::wostringstream strout;
    strout << r1;
    const std::wstring s{strout.str().substr(0, 4)};
    std::wstringstream str;
    str << s;
    str >> r1;
    BOOST_TEST(not str.good(), "stream is not good");
    BOOST_TEST(r1 == r2, "engine not changed");
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

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
  const std::string &name{type_name<T>::name};
  generator_min<R> r_min;
  generator_max<R> r_max;
  using u01xx_min = trng::utility::u01xx_traits<T, 1, generator_min<R>>;
  using u01xx_max = trng::utility::u01xx_traits<T, 1, generator_max<R>>;
  T x_min{}, x_max{};
  x_min = u01xx_min::cc(r_min);
  x_max = u01xx_max::cc(r_max);
  BOOST_TEST((x_min >= 0 && x_max <= 1), "[0, 1] bounds ok for " + name);
  x_min = u01xx_min::co(r_min);
  x_max = u01xx_max::co(r_max);
  BOOST_TEST((x_min >= 0 and x_max < 1), "[0, 1) bounds ok for " + name);
  x_min = u01xx_min::oc(r_min);
  x_max = u01xx_max::oc(r_max);
  BOOST_TEST((x_min > 0 and x_max <= 1), "(0, 1] bounds ok for " + name);
  x_min = u01xx_min::oo(r_min);
  x_max = u01xx_max::oo(r_max);
  BOOST_TEST((x_min > 0 and x_max < 1), "(0, 1) bounds ok for " + name);
}

BOOST_AUTO_TEST_SUITE(test_suite_ranges)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_ranges, R, engines) {
  test_ranges_impl<R, float>();
  test_ranges_impl<R, double>();
  test_ranges_impl<R, long double>();
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_discard)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_discard, R, engines) {
  // two engines with equal state
  R r1, r2;
  const unsigned long long throw_away_1{2ull * 3 * 5 * 7 * 11 * 13 * 17 * 19 + 0x10000000ull};
  const unsigned long long throw_away_2{2ull * 3 * 5 * 7 * 11 * 13 * 17 * 19};
  r1.discard(throw_away_1);
  for (unsigned long long i{0}; i < throw_away_1; ++i)
    r2();
  {
    std::stringstream message;
    message << "engines have equal state after discard(" << throw_away_1 << ")";
    BOOST_TEST(r1 == r2, message.str().c_str());
  }
  r1.discard(throw_away_2);
  for (unsigned long long i{0}; i < throw_away_2; ++i)
    r2();
  {
    std::stringstream message;
    message << "engines have equal state after discard(" << throw_away_2 << ")";
    BOOST_TEST(r1 == r2, message.str().c_str());
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_jump2)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_jump2, R, parallel_engines) {
  // two engines with equal state
  R r1, r2;
  long n{1l};
  // jump ahead in steps of 2^i or step forwards in unit steps
  for (long i{0l}; i < 20l; ++i) {
    r1.jump2(i);
    for (long j{0l}; j < n; ++j)
      r2();
    std::stringstream message;
    message << "engines have equal state after jump2(" << i << ")";
    BOOST_TEST(r1 == r2, message.str().c_str());
    n <<= 1;
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_split)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_split, R, parallel_engines) {
  for (long i{2l}; i <= 20l; ++i) {
    for (long j{0l}; j < i; ++j) {
      // two engines with equal state
      R r1, r2;
      // split or step forwards in j unit steps
      for (long k{0l}; k < j; ++k)
        r1();
      r2.split(i, j);  // split into i streams, pick stream j
      {
        auto v = generate_list(r1, r2, 32, i - 1, 0);
        std::stringstream message;
        message << "engines yield same values after split(" << i << ", " << j << ")";
        BOOST_TEST(std::get<0>(v) == std::get<1>(v), message.str().c_str());
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
