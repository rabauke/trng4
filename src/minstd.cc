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

#include <trng/minstd.hpp>

namespace trng {

  // Parameter and status classes

  // Uniform random number generator concept

  // Equality comparable concept
  bool operator==(const minstd::status_type &S1, 
		  const minstd::status_type &S2) {
    return S1.r==S2.r;
  }

  bool operator!=(const minstd::status_type &S1, 
		  const minstd::status_type &S2) {
    return not (S1==S2);
  }

  // Random number engine concept
  minstd::minstd() : S() { }

  minstd::minstd(unsigned long s) : S() { 
    seed(s);
  }

  void minstd::seed() {
    (*this)=minstd();
  }

  void minstd::seed(unsigned long s) {
    S.r=s%2147483647;
    if (S.r==0)
      S.r=1;
  }

  // Equality comparable concept
  bool operator==(const minstd &R1, const minstd &R2) {
    return R1.S==R2.S;
  }

  bool operator!=(const minstd &R1, const minstd &R2) {
    return not (R1==R2);
  }

  // Other useful methods
  const char * const minstd::name_str="minstd";

  const char * minstd::name() {
    return name_str;
  }

}
