// Copyright (c) 2000-2014, Heiko Bauke
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
#include <fstream>
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

struct parallel_prng { };
struct sequential_prng { };

template<typename T>
struct rng_traits { typedef sequential_prng type; };

template<>
struct rng_traits<trng::lcg64> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::lcg64_shift> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg2> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg3> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg3s> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg4> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg5> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::mrg5s> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn2> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn3> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn3s> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn4> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn5> { typedef parallel_prng type; };

template<>
struct rng_traits<trng::yarn5s> { typedef parallel_prng type; };


template<typename R>
class generator_max {
public:
  typedef typename R::result_type result_type;
#if __cplusplus>=201103L
  static constexpr result_type min() {  return R::min();  }
  static constexpr result_type max() {  return R::max();  }
  result_type operator()() const {
    return max();
  }
#else
  static const result_type min=R::min;
  static const result_type max=R::max;
  result_type operator()() const {
    return max;
  }
#endif
};


template<typename R>
class generator_min {
public:
  typedef typename R::result_type result_type;
#if __cplusplus>=201103L
  static constexpr result_type min() {  return R::min();  }
  static constexpr result_type max() {  return R::max();  }
  result_type operator()() const {
    return min();
  }
#else
  static const result_type min=R::min;
  static const result_type max=R::max;
  result_type operator()() const {
    return min;
  }
#endif
};


template<typename R>
bool test_savestatus_loadstatus(bool genok) {
  // test savestatus / loadstatus 
  R ra, rb;
  for (long i=0l; i<271828l; ++i) {
    ra();
    rb();
  }
  R rc;
  rc=ra;
  for (long i=0l; i<314159l; ++i)
    ra();
  ra=rc;
  bool err=false;
  if (ra()!=rb())
    err=true;
  else
    if (ra()!=rb())
      err=true;
    else
      if (ra()!=rb())
	err=true;
  if (err) {
    std::cout << R::name() 
	      << ": error in savestatus or loadstatus" 
	      << std::endl;
    genok=false;
  }
  return genok;
}


template<typename R>
bool test_status_io(bool genok) {
  // test status i/o
  R ra, rb;
  rb(); rb();
  {
    std::ofstream out("rand.dat");
    out << ra << '\n';
    out.close();
  }
  {
    std::ifstream in("rand.dat");
    in >> rb;
    in.close();
  }
  bool err=false;
  if (ra()!=rb())
    err=true;
  else
    if (ra()!=rb())
      err=true;
    else
      if (ra()!=rb())
	err=true;
  if (err) {
    std::cout << R::name()
	      << ": error in status i/o" << std::endl; 
    genok=false;
  }
  return genok;
}


template<typename R, typename FLOAT>
bool test_ranges_extremes(bool genok) {
  FLOAT x;
  generator_min<R> r_min;
  generator_max<R> r_max;

  x=trng::utility::u01xx_traits<FLOAT, 1, generator_min<R> >::cc(r_min);
  if (x<0) {
    std::cout << R::name()
  	      << ": out of range cc(0)  x = " << x << std::endl; 
    return false;
  }
  x=trng::utility::u01xx_traits<FLOAT, 1, generator_max<R> >::cc(r_max);
  if (x>1) {
    std::cout << R::name()
  	      << ": out of range cc(max)  x = " << x << std::endl; 
    return false;
  }
  
  x=trng::utility::u01xx_traits<FLOAT, 1, generator_min<R> >::co(r_min);
  if (x<0) {
    std::cout << R::name()
  	      << ": out of range co(0)  x = " << x << std::endl; 
    return false;
  }
  x=trng::utility::u01xx_traits<FLOAT, 1, generator_max<R> >::co(r_max);
  if (x>=1) {
    std::cout << R::name()
  	      << ": out of range co(max)  x = " << x << std::endl; 
    return false;
  }

  x=trng::utility::u01xx_traits<FLOAT, 1, generator_min<R> >::oc(r_min);
  if (x<=0) {
    std::cout << R::name()
  	      << ": out of range oc(0)  x = " << x << std::endl; 
    return false;
  }
  x=trng::utility::u01xx_traits<FLOAT, 1, generator_max<R> >::oc(r_max);
  if (x>1) {
    std::cout << R::name()
  	      << ": out of range oc(max)  x = " << x << std::endl; 
    return false;
  }

  x=trng::utility::u01xx_traits<FLOAT, 1, generator_min<R> >::oo(r_min);
  if (x<=0) {
    std::cout << R::name()
  	      << ": out of range oo(0)  x = " << x << std::endl; 
    return false;
  }
  x=trng::utility::u01xx_traits<FLOAT, 1, generator_max<R> >::oo(r_max);
  if (x>=1) {
    std::cout << R::name()
  	      << ": out of range oo(max)  x = " << x << std::endl; 
    return false;
  }

  return genok;
}


template<typename R>
bool test_ranges(bool genok) {
  // test ranges
  genok=test_ranges_extremes<R, float>(genok);
  if (not genok) 
    std::cout << R::name() << ": error for float\n";
  genok=test_ranges_extremes<R, double>(genok);
  if (not genok)
    std::cout << R::name() << ": error for double\n";
  genok=test_ranges_extremes<R, long double>(genok);
  if (not genok)
    std::cout << R::name() << ": error for long double\n";
  return genok;
}


template<typename R>
bool test_jump2(bool genok) {
  // test jump2
  R ra, rb;
  long n=1l;
  for (long i=0l; i<20l; ++i) {
    ra.jump2(i);
    for (long j=0l; j<n; ++j)
      rb();
    bool err=false;
    if (ra()!=rb())
      err=true;
    else
      if (ra()!=rb())
	err=true;
      else
	if (ra()!=rb())
	  err=true;
    if (err) {
      std::cout << ra.name() << ": error in jump2" << std::endl;
      genok=false;
      break;
      }
    n<<=1;
  }
  return genok;
}


template<typename R>
bool test_split(bool genok) {
  // test split
  for (long i=2l; i<=20l; ++i) {
    for (long j=0l; j<i; ++j) {
      R ra, rb;
      for (long k=0l; k<j; ++k)
	ra();
      rb.split(i, j);
      bool err=false;
      if (ra()!=rb())
	err=true;
      else {
	for (long k=0l; k<i-1l; ++k)
	  ra();
	if (ra()!=rb())
	  err=true;
	else {
	  for (long k=0l; k<i-1l; ++k)
	    ra();
	  if (ra()!=rb())
	    err=true;
	}
      }
      if (err) {
	std::cout << ra.name() << ": error in split" << std::endl;
	genok=false;
	break;
      }
    }
  }
  return genok;
}


template<typename R>
void plausibility_main_impl(sequential_prng) {
  bool genok(true);
  std::cout << "testing sequential PRNG " << R::name() << std::endl;
  genok=test_savestatus_loadstatus<R>(genok);
  genok=test_status_io<R>(genok);
  genok=test_ranges<R>(genok);
  std::cout << R::name() << ": test " 
	    << (genok ? "passed" : "failed") << std::endl;
}


template<typename R>
void plausibility_main_impl(parallel_prng) {
  bool genok(true);
  std::cout << "testing parallel PRNG " << R::name() << std::endl;
  genok=test_savestatus_loadstatus<R>(genok);
  genok=test_status_io<R>(genok);
  genok=test_ranges<R>(genok);
  genok=test_jump2<R>(genok);
  genok=test_split<R>(genok);
  std::cout << R::name() << ": test " 
	    << (genok ? "passed" : "failed") << std::endl;
}


template<typename R>
void plausibility_main() {
  plausibility_main_impl<R>(typename rng_traits<R>::type());
}


int main() {
  try {
    { plausibility_main<trng::lcg64>(); }
    { plausibility_main<trng::lcg64_shift>(); }
    { plausibility_main<trng::mrg2>(); }
    { plausibility_main<trng::mrg3>(); }
    { plausibility_main<trng::mrg3s>(); }
    { plausibility_main<trng::mrg4>(); }
    { plausibility_main<trng::mrg5>(); }
    { plausibility_main<trng::mrg5s>(); }
    { plausibility_main<trng::yarn2>(); }
    { plausibility_main<trng::yarn3>(); }
    { plausibility_main<trng::yarn3s>(); }
    { plausibility_main<trng::yarn4>(); }
    { plausibility_main<trng::yarn5>(); }
    { plausibility_main<trng::yarn5s>(); }
    { plausibility_main<trng::lagfib2xor_521_32>(); }
    { plausibility_main<trng::lagfib4xor_521_32>(); }
    { plausibility_main<trng::lagfib2plus_521_32>(); }
    { plausibility_main<trng::lagfib4plus_521_32>(); }
    { plausibility_main<trng::mt19937>(); }
    { plausibility_main<trng::mt19937_64>(); }
  } 
  catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
  }
  catch (...) {
    std::cerr << "something went wrong" << std::endl;
  }
  return EXIT_SUCCESS;
}


