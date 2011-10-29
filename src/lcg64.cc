// Copyright (C) 2000-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
//  
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License in
// version 2 as published by the Free Software Foundation.
//  
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
//  

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
    return !(P1==P2);
  }
  
  // Equality comparable concept
  bool operator==(const lcg64::status_type &S1, 
		  const lcg64::status_type &S2) {
    return S1.r==S2.r;
  }

  bool operator!=(const lcg64::status_type &S1, 
		  const lcg64::status_type &S2) {
    return !(S1==S2);
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
    return !(R1==R2);
  }

  // Parallel random number generator concept
  void lcg64::split(unsigned int s, unsigned int n) {
    if (s<1 or n>=s)
      throw std::invalid_argument("invalid argument for trng::lcg64::split");
    if (s>1) {
      lcg64::result_type t1(1ull), t2(0ull);
      for (unsigned int i(0); i<=n; ++i)
	step();
      for (unsigned int i(0); i<s; ++i) {
	t2+=t1;
	t1*=P.a;
      }
      P.b*=t2;
      P.a=t1;
      backward();
    }
  }

  void lcg64::jump2(unsigned int s) {
    lcg64::result_type t1(P.a);
    for (unsigned int i(0); i<s; ++i)
      t1*=t1;
    lcg64::result_type t2(1ull), t3(P.a);
    while (s>0l) {
      t2*=(1ull+t3);
      t3*=t3;
      --s;
    }
    S.r=S.r*t1+t2*P.b;
  }

  void lcg64::jump(unsigned long long s) {
    unsigned int i(0);
    while (s>0) {
      if (s%2u==1u)
	jump2(i);
      ++i;
      s>>=1;
    }
  }
  
  // Other usefull methods
  const char * const lcg64::name_str="lcg64";
  
  const char * lcg64::name() {
    return name_str;
  }
  
  void lcg64::backward() {
    for (unsigned int i(0); i<64; ++i)
      jump2(i);
  }

}

