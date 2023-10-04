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

#include <cstdlib>
#include <iostream>
#include <exception>
#include <string>
#include <sstream>
#include <chrono>
#include <random>
#include <trng/generate_canonical.hpp>
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
#include <trng/mt19937_64.hpp>
#include <trng/mt19937.hpp>
#include <trng/lagfib2xor.hpp>
#include <trng/lagfib2plus.hpp>
#include <trng/lagfib4xor.hpp>
#include <trng/lagfib4plus.hpp>
#include <trng/xoshiro256plus.hpp>

#if defined TRNG_HAVE_BOOST
#include <boost/random/linear_congruential.hpp>
#include <boost/random/additive_combine.hpp>
#include <boost/random/inversive_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/shuffle_output.hpp>
#endif

template<typename T>
std::string to_string(const T &x) {
  std::ostringstream temp;
  temp << x;
  return temp.str();
}

class timer {
private:
  const double _resolution;
  std::chrono::time_point<std::chrono::high_resolution_clock> _t;

public:
  void reset() { _t = std::chrono::high_resolution_clock::now(); }
  double time() const {
    const auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(now - _t).count() * 1e-6;
  }
  double resolution() const { return _resolution; }
  timer()
      : _resolution([]() { return 1e-6; }()), _t(std::chrono::high_resolution_clock::now()) {}
};

template<typename R>
typename R::result_type time_plain(R &r) {
  typename R::result_type s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += r();
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R>
double time_cc(R &r) {
  double s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += trng::utility::uniformcc<double>(r);
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R>
double time_oc(R &r) {
  double s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += trng::utility::uniformoc<double>(r);
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R>
double time_co(R &r) {
  double s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += trng::utility::uniformco<double>(r);
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R>
double time_oo(R &r) {
  double s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += trng::utility::uniformoo<double>(r);
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R, typename FLOAT>
FLOAT time_canonical(R &r) {
  FLOAT s{0};
  timer T;
  long max(1l << 24);
  for (long i{0}; i < max; ++i)
    s += trng::generate_canonical<FLOAT>(r);
  std::string res(to_string(1e-6 * max / T.time()));
  while (res.length() < 10)
    res += ' ';
  std::cout << res;
  return s;
}

template<typename R>
void time_main(R &r, std::string name) {
  std::stringstream s;  // write data to stream to prevent that code gets optimized away
  while (name.length() < 32)
    name += ' ';
  std::cout << name;
  const auto res_plain = time_plain(r);
  s << res_plain;
  const auto res_cc = time_cc(r);
  s << res_cc;
  const auto res_co = time_co(r);
  s << res_co;
  const auto res_oc = time_oc(r);
  s << res_oc;
  const auto res_oo = time_oo(r);
  s << res_oo;
  const auto res_canonical = time_canonical<R, double>(r);
  s << res_canonical;
  std::cout << std::endl;
}

template<typename R>
void time_boost(R &r, std::string name) {
  std::stringstream s;  // write data to stream to prevent that code gets optimized away
  while (name.length() < 32)
    name += ' ';
  std::cout << name;
  const auto res_plain = time_plain(r);
  s << res_plain;
  std::cout << std::endl;
}

int main() {
  std::cout << "                                            10^6 random numbers per second\n"
            << "generator                       [min,max] [0,1]     [0,1)     (0,1]     (0,1)  "
               "   canonical\n"
            << "==============================================================================="
               "==============\n";
  std::cout.flush();
  try {
    {
      trng::lcg64 r;
      time_main(r, "trng::lcg64");
    }
    {
      trng::lcg64_shift r;
      time_main(r, "trng::lcg64_shift");
    }
    {
      trng::lcg64_count_shift r;
      time_main(r, "trng::lcg64_count_shift");
    }
    {
      trng::mrg2 r;
      time_main(r, "trng::mrg2");
    }
    {
      trng::mrg3 r;
      time_main(r, "trng::mrg3");
    }
    {
      trng::mrg3s r;
      time_main(r, "trng::mrg3s");
    }
    {
      trng::mrg4 r;
      time_main(r, "trng::mrg4");
    }
    {
      trng::mrg5 r;
      time_main(r, "trng::mrg5");
    }
    {
      trng::mrg5s r;
      time_main(r, "trng::mrg5s");
    }
    {
      trng::yarn2 r;
      time_main(r, "trng::yarn2");
    }
    {
      trng::yarn3 r;
      time_main(r, "trng::yarn3");
    }
    {
      trng::yarn3s r;
      time_main(r, "trng::yarn3s");
    }
    {
      trng::yarn4 r;
      time_main(r, "trng::yarn4");
    }
    {
      trng::yarn5 r;
      time_main(r, "trng::yarn5");
    }
    {
      trng::yarn5s r;
      time_main(r, "trng::yarn5s");
    }
    {
      trng::mt19937 r;
      time_main(r, "trng::mt19937");
    }
    {
      trng::mt19937_64 r;
      time_main(r, "trng::mt19937_64");
    }
    {
      trng::lagfib2xor_19937_64 r;
      time_main(r, "trng::lagfib2xor_19937_64");
    }
    {
      trng::lagfib4xor_19937_64 r;
      time_main(r, "trng::lagfib4xor_19937_64");
    }
    {
      trng::lagfib2plus_19937_64 r;
      time_main(r, "trng::lagfib2plus_19937_64");
    }
    {
      trng::lagfib4plus_19937_64 r;
      time_main(r, "trng::lagfib4plus_19937_64");
    }
    {
      trng::xoshiro256plus r;
      time_main(r, "trng::xoshiro256plus");
    }
    {
      std::minstd_rand0 r;
      time_main(r, "std::minstd_rand0");
    }
    {
      std::minstd_rand r;
      time_main(r, "std::minstd_rand");
    }
    {
      std::mt19937 r;
      time_main(r, "std::mt19937");
    }
    {
      std::mt19937_64 r;
      time_main(r, "std::mt19937_64");
    }
    {
      std::ranlux24_base r;
      time_main(r, "std::ranlux24_base");
    }
    {
      std::ranlux48_base r;
      time_main(r, "std::ranlux48_base");
    }
    {
      std::ranlux24 r;
      time_main(r, "std::ranlux24");
    }
    {
      std::ranlux48 r;
      time_main(r, "std::ranlux48");
    }
    {
      std::knuth_b r;
      time_main(r, "std::knuth_b");
    }
#if defined TRNG_HAVE_BOOST
    {
      boost::minstd_rand r;
      time_boost(r, "boost::minstd_rand");
    }
    {
      boost::ecuyer1988 r;
      time_boost(r, "boost::ecuyer1988");
    }
    {
      boost::kreutzer1986 r;
      time_boost(r, "boost::kreutzer1986");
    }
    {
      boost::hellekalek1995 r;
      time_boost(r, "boost::hellekalek1995");
    }
    {
      boost::mt11213b r;
      time_boost(r, "boost::mt11213b");
    }
    {
      boost::mt19937 r;
      time_boost(r, "boost::mt19937");
    }
    {
      boost::lagged_fibonacci607 r;
      time_boost(r, "boost::lagged_fibonacci607");
    }
    {
      boost::lagged_fibonacci1279 r;
      time_boost(r, "boost::lagged_fibonacci1279");
    }
    {
      boost::lagged_fibonacci2281 r;
      time_boost(r, "boost::lagged_fibonacci2281");
    }
    {
      boost::lagged_fibonacci3217 r;
      time_boost(r, "boost::lagged_fibonacci3217");
    }
    {
      boost::lagged_fibonacci4423 r;
      time_boost(r, "boost::lagged_fibonacci4423");
    }
    {
      boost::lagged_fibonacci9689 r;
      time_boost(r, "boost::lagged_fibonacci9689");
    }
    {
      boost::lagged_fibonacci19937 r;
      time_boost(r, "boost::lagged_fibonacci19937");
    }
    {
      boost::lagged_fibonacci23209 r;
      time_boost(r, "boost::lagged_fibonacci23209");
    }
    {
      boost::lagged_fibonacci44497 r;
      time_boost(r, "boost::lagged_fibonacci44497");
    }
#endif
  } catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
  } catch (...) {
    std::cerr << "something went wrong" << std::endl;
  }
  return EXIT_SUCCESS;
}
