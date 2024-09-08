// Copyright (c) 2000-2024, Heiko Bauke
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

#include "count128_lcg_shift.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // Equality comparable concept
  bool operator==(const count128_lcg_shift::parameter_type &P1,
                  const count128_lcg_shift::parameter_type &P2) {
    return P1.a == P2.a and P1.b == P2.b;
  }

  bool operator!=(const count128_lcg_shift::parameter_type &P1,
                  const count128_lcg_shift::parameter_type &P2) {
    return not(P1 == P2);
  }

  // Equality comparable concept
  bool operator==(const count128_lcg_shift::status_type &S1, const count128_lcg_shift::status_type &S2) {
    return S1.r == S2.r;
  }

  bool operator!=(const count128_lcg_shift::status_type &S1, const count128_lcg_shift::status_type &S2) {
    return not(S1 == S2);
  }

  // 128-bit increment equals the prime-number 337796325545380861827125810166389624843
  const count128_lcg_shift::parameter_type count128_lcg_shift::Default =
      parameter_type(uint128{0xfe2134b266a61770, 0x32095479a8f5500b}, 18145460002477866997u, 1u);
  const count128_lcg_shift::parameter_type count128_lcg_shift::LEcuyer1 =
      parameter_type(uint128{0xfe2134b266a61770, 0x32095479a8f5500b}, 2862933555777941757u, 1u);
  const count128_lcg_shift::parameter_type count128_lcg_shift::LEcuyer2 =
      parameter_type(uint128{0xfe2134b266a61770, 0x32095479a8f5500b}, 3202034522624059733u, 1u);
  const count128_lcg_shift::parameter_type count128_lcg_shift::LEcuyer3 =
      parameter_type(uint128{0xfe2134b266a61770, 0x32095479a8f5500b}, 3935559000370003845u, 1u);

  // Random number engine concept
  count128_lcg_shift::count128_lcg_shift(count128_lcg_shift::parameter_type P) : P{P} {}

  count128_lcg_shift::count128_lcg_shift(unsigned long s, count128_lcg_shift::parameter_type P) : P{P} { seed(s);
  }

  count128_lcg_shift::count128_lcg_shift(unsigned long long s,
                                         count128_lcg_shift::parameter_type P) : P{P} {
    seed(s);
  }

  void count128_lcg_shift::seed() { (*this) = count128_lcg_shift(); }

  void count128_lcg_shift::seed(unsigned long s) {
    S.r = uint128{static_cast<count128_lcg_shift::result_type>(s), static_cast<count128_lcg_shift::result_type>(s)};
  }

  void count128_lcg_shift::seed(unsigned long long s) {
    S.r = uint128{static_cast<count128_lcg_shift::result_type>(s), static_cast<count128_lcg_shift::result_type>(s)};
  }

  // Equality comparable concept
  bool operator==(const count128_lcg_shift &R1, const count128_lcg_shift &R2) {
    return R1.P == R2.P and R1.S == R2.S;
  }

  bool operator!=(const count128_lcg_shift &R1, const count128_lcg_shift &R2) { return not(R1 == R2); }

  // Parallel random number generator concept

  // Other useful methods
  const char *const count128_lcg_shift::name_str = "count128_lcg_shift";

  const char *count128_lcg_shift::name() { return name_str; }

}  // namespace trng
