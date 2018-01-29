// Copyright (c) 2000-2018, Heiko Bauke
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

#include <trng/mt19937_64.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const mt19937_64::parameter_type &, 
		  const mt19937_64::parameter_type &) {
    return true;
  }

  bool operator!=(const mt19937_64::parameter_type &, 
		  const mt19937_64::parameter_type &) {
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
    return not (S1==S2);
  }
  
  // Random number engine concept
  mt19937_64::mt19937_64() :
    P(), S() { 
    seed(5489u); 
  }

  mt19937_64::mt19937_64(unsigned long s) :
    P(), S() { 
    seed(s);
  }
    
  void mt19937_64::seed() {
    (*this)=mt19937_64();
  }
 
  void mt19937_64::seed(unsigned long s) {
    S.mt[0]=s;
    for (S.mti=1; S.mti<mt19937_64::status_type::N; ++S.mti) 
      S.mt[S.mti]=(6364136223846793005u * (S.mt[S.mti-1] ^ (S.mt[S.mti-1] >> 62)) + S.mti);
  }
  
  // Equality comparable concept
  bool operator==(const mt19937_64 &R1, const mt19937_64 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const mt19937_64 &R1, const mt19937_64 &R2) {
    return not (R1==R2);
  }

  // Other useful methods
  const char * const mt19937_64::name_str="mt19937_64";
  
  const char * mt19937_64::name() {
    return name_str;
  }
  
}

