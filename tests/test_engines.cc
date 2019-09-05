#define BOOST_TEST_MODULE test engines
#define BOOST_TEST_DYN_LINK

#include <vector>
#include <sstream>
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
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

typedef boost::mpl::list<trng::lcg64, trng::lcg64_shift,  //
                         trng::mrg2, trng::mrg3, trng::mrg3s, trng::mrg4, trng::mrg5,
                         trng::mrg5s,  //
                         trng::yarn2, trng::yarn3, trng::yarn3s, trng::yarn4, trng::yarn5,
                         trng::yarn5s,                                        //
                         trng::lagfib2xor_521_64, trng::lagfib4xor_521_32,    //
                         trng::lagfib2plus_521_32, trng::lagfib4plus_521_64,  //
                         trng::mt19937, trng::mt19937_64>
    engines;


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
generate_list(R &r1, R &r2, const long N) {
  std::vector<typename R::result_type> v1, v2;
  if (N > 0) {
    v1.reserve(N);
    v2.reserve(N);
  }
  for (long i{0}; i < N; ++i) {
    v1.push_back(r1());
    v2.push_back(r2());
  }
  return {v1, v2};
}


BOOST_AUTO_TEST_SUITE(test_suite_engines)

BOOST_AUTO_TEST_SUITE(test_suite_advance)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_advance, R, engines) {
  // two generators with equal state
  R r1, r2;
  advance_engine(r1, 271828l);  // advance generator ra
  BOOST_TEST(r1 != r2, "engines have different state");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_suite_restore)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_restore, R, engines) {
  // three generators with equal state
  R r1, r2, r3;
  advance_engine(r1, r2, 271828l);  // advance generators ra and rb
  r3 = r1;                          // backup ra in rc
  // advance ra such that ra and rb have different states
  advance_engine(r1, 314159l);
  r1 = r3;  // restore ra
  BOOST_TEST(r1 == r2, "engine state restored successfully");
  // when state has been restored correctly, ra and rb must yield the same values
  auto v = generate_list(r1, r2, 32);
  BOOST_TEST(std::get<0>(v) == std::get<1>(v), "engines yield same values after restore");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_suite_status_io)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_status_io, R, engines) {
  // two generators with equal state
  R r1, r2;
  advance_engine(r1, 271828l);  // advance generator ra
  std::stringstream str;
  str << r1 << '\n';
  str >> r2;
  BOOST_TEST(r1 == r2, "engine state restored from stream successfully");
  // when state has been restored correctly, ra and rb must yield the same values
  auto v = generate_list(r1, r2, 32);
  BOOST_TEST(std::get<0>(v) == std::get<1>(v),
             "engines yield same values after restore from stream");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
