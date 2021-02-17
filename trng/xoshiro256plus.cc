// Copyright (c) 2000-2021, Heiko Bauke
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

#include "xoshiro256plus.hpp"
#include "minstd.hpp"

namespace trng {

  // Uniform random number generator concept

  // Status class

  // Equality comparable concept
  bool operator==(const xoshiro256plus::status_type &S1,
                  const xoshiro256plus::status_type &S2) {
    return S1.r[0] == S2.r[0] and S1.r[1] == S2.r[1] and S1.r[2] == S2.r[2] and
           S1.r[3] == S2.r[3];
  }

  bool operator!=(const xoshiro256plus::status_type &S1,
                  const xoshiro256plus::status_type &S2) {
    return not(S1 == S2);
  }

  // Random number engine concept
  xoshiro256plus::xoshiro256plus() = default;

  xoshiro256plus::xoshiro256plus(unsigned long s) { seed(s); }

  xoshiro256plus::xoshiro256plus(result_type s0, result_type s1, result_type s2,
                                 result_type s3) {
    seed(s0, s1, s2, s3);
  }

  void xoshiro256plus::seed() { (*this) = xoshiro256plus(); }

  void xoshiro256plus::seed(unsigned long s) {
    minstd R(s);
    seed(R);
  }

  void xoshiro256plus::seed(result_type s0, result_type s1, result_type s2, result_type s3) {
    S.r[0] = static_cast<result_type>(s0);
    S.r[1] = static_cast<result_type>(s1);
    S.r[2] = static_cast<result_type>(s2);
    S.r[3] = static_cast<result_type>(s3);
    if (S.r[0] == 0 and S.r[1] == 0 and S.r[2] == 0 and S.r[3] == 0)
      S.r[0] |= result_type(1) << 63u;
  }

  // Equality comparable concept
  bool operator==(const xoshiro256plus &R1, const xoshiro256plus &R2) { return R1.S == R2.S; }

  bool operator!=(const xoshiro256plus &R1, const xoshiro256plus &R2) { return not(R1 == R2); }

  // Parallel random number generator concept

  // Other useful methods
  const char *const xoshiro256plus::name_str = "xoshiro256plus";

  const char *xoshiro256plus::name() { return name_str; }

}  // namespace trng
