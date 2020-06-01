// Copyright (c) 2000-2020, Heiko Bauke
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

#include "mrg2.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes
  const mrg2::parameter_type mrg2::LEcuyer1 = parameter_type(1498809829, 1160990996);
  const mrg2::parameter_type mrg2::LEcuyer2 = parameter_type(46325, 1084587);

  // Random number engine concept
  mrg2::mrg2(mrg2::parameter_type P) : P{P} {}

  mrg2::mrg2(unsigned long s, mrg2::parameter_type P) : P{P} { seed(s); }

  void mrg2::seed() { (*this) = mrg2(); }

  void mrg2::seed(unsigned long s) {
    int64_t t(s);
    t %= modulus;
    if (t < 0)
      t += modulus;
    S.r[0] = static_cast<result_type>(t);
    S.r[1] = 1;
  }

  void mrg2::seed(mrg2::result_type s1, mrg2::result_type s2) {
    S.r[0] = s1 % modulus;
    if (S.r[0] < 0)
      S.r[0] += modulus;
    S.r[1] = s2 % modulus;
    if (S.r[1] < 0)
      S.r[1] += modulus;
  }

  // Equality comparable concept
  bool operator==(const mrg2 &R1, const mrg2 &R2) { return R1.P == R2.P and R1.S == R2.S; }

  bool operator!=(const mrg2 &R1, const mrg2 &R2) { return not(R1 == R2); }

  // other useful methods
  const char *const mrg2::name_str = "mrg2";

  const char *mrg2::name() { return name_str; }

}  // namespace trng
