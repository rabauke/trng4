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

#include <trng/mt19937_64.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const mt19937_64::parameter_type &P1, 
		  const mt19937_64::parameter_type &P2) {
    return true;
  }

  bool operator!=(const mt19937_64::parameter_type &P1, 
		  const mt19937_64::parameter_type &P2) {
    return false;
  }
  
  // Equality comparable concept
  bool operator==(const mt19937_64::status_type &S1, 
		  const mt19937_64::status_type &S2) {
    for (int i=0; i<mt19937_64::status_type::N; ++i)
      if (S1.mt[i]!=S2.mt[i])
	return false;
    return true;
  }

  bool operator!=(const mt19937_64::status_type &S1, 
		  const mt19937_64::status_type &S2) {
    return !(S1==S2);
  }
  
  // Random number engine concept
  mt19937_64::mt19937_64() :
    P(), S() { 
    seed(5489ULL); 
  }

  mt19937_64::mt19937_64(unsigned long s) :
    P(), S() { 
    seed(s);
  }
    
  void mt19937_64::seed() {
    (*this)=mt19937_64();
  }
 
  void mt19937_64::seed(unsigned long s) {
    seed(static_cast<mt19937_64::result_type>(s));
  }
  
  void mt19937_64::seed(mt19937_64::result_type s) {
    S.mt[0]=s;
    for (S.mti=1; S.mti<mt19937_64::status_type::N; ++S.mti) 
      S.mt[S.mti]=(6364136223846793005ull * (S.mt[S.mti-1] ^ (S.mt[S.mti-1] >> 62)) + S.mti);
  }
  
  // Equality comparable concept
  bool operator==(const mt19937_64 &R1, const mt19937_64 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const mt19937_64 &R1, const mt19937_64 &R2) {
    return !(R1==R2);
  }

  // Other usefull methods
  const char * const mt19937_64::name_str="mt19937_64";
  
  const char * mt19937_64::name() {
    return name_str;
  }
  
}

