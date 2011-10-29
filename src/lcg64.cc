// Copyright (c) 2000-2010, Heiko Bauke
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

#include <trng/lcg64.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const lcg64::parameter_type &P1, 
		  const lcg64::parameter_type &P2) {
    return P1.a==P2.a and P1.b==P2.b;
  }

  bool operator!=(const lcg64::parameter_type &P1, 
		  const lcg64::parameter_type &P2) {
    return not (P1==P2);
  }
  
  // Equality comparable concept
  bool operator==(const lcg64::status_type &S1, 
		  const lcg64::status_type &S2) {
    return S1.r==S2.r;
  }

  bool operator!=(const lcg64::status_type &S1, 
		  const lcg64::status_type &S2) {
    return not (S1==S2);
  }
  
  // compute floor(log_2(x))
  unsigned int lcg64::log2_floor(lcg64::result_type x) {
    unsigned int y(0);
    while (x>0) {
      x>>=1;
      ++y;
    };
    --y;
    return y;
  }

  // compute pow(x, n)
  lcg64::result_type lcg64::pow(lcg64::result_type x, lcg64::result_type n) {
    lcg64::result_type result=1;
    while (n>0) {
      if ((n&1)>0)
	result=result*x;
      x=x*x;
      n>>=1;
    }
    return result;
  }

  // compute prod(1+a^(2^i), i=0..l-1)
  lcg64::result_type lcg64::g(unsigned int l, lcg64::result_type a) {
    lcg64::result_type p=a, res=1;
    for (unsigned i=0; i<l; ++i) {
      res*=1+p;
      p*=p;
    }
#if ULONG_LONG_MAX>18446744073709551615ull
    return res & 0xfffffffffffffffful;
#else
    return res;
#endif
  }
  
  // compute sum(a^i, i=0..s-1)
  lcg64::result_type lcg64::f(lcg64::result_type s, lcg64::result_type a) {
    if (s==0)
      return 0;
    unsigned int l=log2_floor(s);
    lcg64::result_type p=a;
    for (unsigned i=0; i<l; ++i)
      p*=p;
#if ULONG_LONG_MAX>18446744073709551615ull
    return ( g(l, a)+p*f(s-(1ull<<l), a) ) & 0xfffffffffffffffful;
#else
    return g(l, a)+p*f(s-(1ull<<l), a);
#endif
  }

  const lcg64::parameter_type
  lcg64::Default=parameter_type(18145460002477866997ull, 1ull);
  const lcg64::parameter_type 
  lcg64::LEcuyer1=parameter_type(2862933555777941757ull, 1ull);
  const lcg64::parameter_type 
  lcg64::LEcuyer2=parameter_type(3202034522624059733ull, 1ull);
  const lcg64::parameter_type 
  lcg64::LEcuyer3=parameter_type(3935559000370003845ull, 1ull);

  // Random number engine concept
  lcg64::lcg64(lcg64::parameter_type P) :
    P(P), S() { }

  lcg64::lcg64(unsigned long s, lcg64::parameter_type P) :
    P(P), S() { 
    seed(s);
  }
    
  void lcg64::seed() {
    (*this)=lcg64();
  }
 
  void lcg64::seed(unsigned long s) {
    S.r=s;
  }
  
  void lcg64::seed(lcg64::result_type s) {
    S.r=s;
#if ULONG_LONG_MAX>18446744073709551615ull
    S.r&=0xfffffffffffffffful;
#endif
  }
  
  // Equality comparable concept
  bool operator==(const lcg64 &R1, const lcg64 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const lcg64 &R1, const lcg64 &R2) {
    return not (R1==R2);
  }

  // Parallel random number generator concept
  void lcg64::split(unsigned int s, unsigned int n) {
    if (s<1 or n>=s)
      throw std::invalid_argument("invalid argument for trng::lcg64::split");
    if (s>1) {
      jump(n+1);
      P.b*=f(s, P.a);
      P.a=pow(P.a, s);
      backward();
    }
  }

  void lcg64::jump2(unsigned int s) {
    S.r=S.r*pow(P.a, 1ull<<s)+f(1ull<<s, P.a)*P.b;
#if ULONG_LONG_MAX>18446744073709551615ull
    S.r&=0xfffffffffffffffful;
#endif
  }

  void lcg64::jump(unsigned long long s) {
    if (s<16) {
      for (unsigned int i(0); i<s; ++i) 
	step();
    } else {
      unsigned int i(0);
      while (s>0) {
	if (s%2==1)
	  jump2(i);
	++i;
	s>>=1;
      }
    }
  }
  
  // Other useful methods
  const char * const lcg64::name_str="lcg64";
  
  const char * lcg64::name() {
    return name_str;
  }
  
  void lcg64::backward() {
    for (unsigned int i(0); i<64; ++i)
      jump2(i);
  }

}

