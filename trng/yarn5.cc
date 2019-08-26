// Copyright (c) 2000-2019, Heiko Bauke
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

#include "yarn5.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // EqualityComparable concept
  bool operator==(const yarn5::parameter_type &P1, const yarn5::parameter_type &P2) {
    return P1.a1 == P2.a1 and P1.a2 == P2.a2 and P1.a3 == P2.a3 and P1.a4 == P2.a4 and
           P1.a5 == P2.a5;
  }

  bool operator!=(const yarn5::parameter_type &P1, const yarn5::parameter_type &P2) {
    return not(P1 == P2);
  }

  // Equality comparable concept
  bool operator==(const yarn5::status_type &S1, const yarn5::status_type &S2) {
    return S1.r1 == S2.r1 and S1.r2 == S2.r2 and S1.r3 == S2.r3 and S1.r4 == S2.r4 and
           S1.r5 == S2.r5;
  }

  bool operator!=(const yarn5::status_type &S1, const yarn5::status_type &S2) {
    return not(S1 == S2);
  }

  const yarn5::parameter_type yarn5::LEcuyer1 = parameter_type(107374182, 0, 0, 0, 104480);

  // Random number engine concept
  yarn5::yarn5(const parameter_type &P) : P(P), S() {}

  yarn5::yarn5(unsigned long s, const parameter_type &P) : P(P), S() { seed(s); }

  void yarn5::seed() { (*this) = yarn5(); }

  void yarn5::seed(unsigned long s) {
    int64_t t = s;
    t %= modulus;
    if (t < 0)
      t += modulus;
    S.r1 = static_cast<result_type>(t);
    S.r2 = 1;
    S.r3 = 1;
    S.r4 = 1;
    S.r5 = 1;
  }

  void yarn5::seed(yarn5::result_type s1, yarn5::result_type s2, yarn5::result_type s3,
                   yarn5::result_type s4, yarn5::result_type s5) {
    S.r1 = s1 % modulus;
    if (S.r1 < 0)
      S.r1 += modulus;
    S.r2 = s2 % modulus;
    if (S.r2 < 0)
      S.r2 += modulus;
    S.r3 = s3 % modulus;
    if (S.r3 < 0)
      S.r3 += modulus;
    S.r4 = s4 % modulus;
    if (S.r4 < 0)
      S.r4 += modulus;
    S.r5 = s5 % modulus;
    if (S.r5 < 0)
      S.r5 += modulus;
  }

  // Equality comparable concept
  bool operator==(const yarn5 &R1, const yarn5 &R2) { return R1.P == R2.P and R1.S == R2.S; }

  bool operator!=(const yarn5 &R1, const yarn5 &R2) { return not(R1 == R2); }

  // other useful methods
  const char *const yarn5::name_str = "yarn5";

  const char *yarn5::name() { return name_str; }

  int_math::power<yarn5::modulus, yarn5::gen> yarn5::parameter_type::g;

}  // namespace trng
