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

#include <trng/mt19937.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const mt19937::parameter_type &P1, 
		  const mt19937::parameter_type &P2) {
    return true;
  }

  bool operator!=(const mt19937::parameter_type &P1, 
		  const mt19937::parameter_type &P2) {
    return false;
  }
  
  // Equality comparable concept
  bool operator==(const mt19937::status_type &S1, 
		  const mt19937::status_type &S2) {
    for (int i=0; i<mt19937::status_type::N; ++i)
      if (S1.mt[i]!=S2.mt[i])
	return false;
    return true;
  }

  bool operator!=(const mt19937::status_type &S1, 
		  const mt19937::status_type &S2) {
    return !(S1==S2);
  }
  
  // Random number engine concept
  mt19937::mt19937() :
    P(), S() { 
    seed(5489ULL); 
  }

  mt19937::mt19937(unsigned long s) :
    P(), S() { 
    seed(s);
  }
    
  void mt19937::seed() {
    (*this)=mt19937();
  }
 
  void mt19937::seed(mt19937::result_type s) {
    S.mt[0]=s & 0xffffffffUL;
    for (S.mti=1; S.mti<mt19937::N; ++S.mti) 
      S.mt[S.mti]=(1812433253UL * (S.mt[S.mti-1] ^ (S.mt[S.mti-1] >> 30)) + S.mti) & 0xffffffffUL; 
  }
  
  // Equality comparable concept
  bool operator==(const mt19937 &R1, const mt19937 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const mt19937 &R1, const mt19937 &R2) {
    return !(R1==R2);
  }

  // Other usefull methods
  const char * const mt19937::name_str="mt19937";
  
  const char * mt19937::name() {
    return name_str;
  }
  
}

