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

#include "lcg64.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const lcg64::parameter_type &P1, const lcg64::parameter_type &P2) {
    return P1.a == P2.a and P1.b == P2.b;
  }

  bool operator!=(const lcg64::parameter_type &P1, const lcg64::parameter_type &P2) {
    return not(P1 == P2);
  }

  // Equality comparable concept
  bool operator==(const lcg64::status_type &S1, const lcg64::status_type &S2) {
    return S1.r == S2.r;
  }

  bool operator!=(const lcg64::status_type &S1, const lcg64::status_type &S2) {
    return not(S1 == S2);
  }

  const lcg64::parameter_type lcg64::Default = parameter_type(18145460002477866997u, 1u);
  const lcg64::parameter_type lcg64::LEcuyer1 = parameter_type(2862933555777941757u, 1u);
  const lcg64::parameter_type lcg64::LEcuyer2 = parameter_type(3202034522624059733u, 1u);
  const lcg64::parameter_type lcg64::LEcuyer3 = parameter_type(3935559000370003845u, 1u);

  // Random number engine concept
  lcg64::lcg64(lcg64::parameter_type P) : P{P} {}

  lcg64::lcg64(unsigned long s, lcg64::parameter_type P) : P{P} { seed(s); }

  lcg64::lcg64(unsigned long long s, lcg64::parameter_type P) : P{P} { seed(s); }

  void lcg64::seed() { (*this) = lcg64(); }

  void lcg64::seed(unsigned long s) { S.r = static_cast<result_type>(s); }

  void lcg64::seed(unsigned long long s) { S.r = static_cast<result_type>(s); }

  // Equality comparable concept
  bool operator==(const lcg64 &R1, const lcg64 &R2) { return R1.P == R2.P and R1.S == R2.S; }

  bool operator!=(const lcg64 &R1, const lcg64 &R2) { return not(R1 == R2); }

  // Parallel random number generator concept

  // Other useful methods
  const char *const lcg64::name_str = "lcg64";

  const char *lcg64::name() { return name_str; }

}  // namespace trng
