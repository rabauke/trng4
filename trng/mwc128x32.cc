// Copyright (c) 2000-2026, Heiko Bauke
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

#include "mwc128x32.hpp"
#include "minstd.hpp"

namespace trng {

  // Uniform random number generator concept

  // Status class

  // Equality comparable concept
  bool operator==(const mwc128x32::status_type &S1, const mwc128x32::status_type &S2) {
    return S1.x == S2.x and S1.y == S2.y and S1.z == S2.z and S1.c == S2.c;
  }

  bool operator!=(const mwc128x32::status_type &S1, const mwc128x32::status_type &S2) {
    return not(S1 == S2);
  }

  // Random number engine concept
  mwc128x32::mwc128x32() = default;

  mwc128x32::mwc128x32(unsigned long s) {
    seed(s);
  }

  mwc128x32::mwc128x32(result_type x, result_type y, result_type z, result_type c) {
    seed(x, y, z, c);
  }

  void mwc128x32::seed() {
    (*this) = mwc128x32();
  }

  void mwc128x32::seed(unsigned long s) {
    minstd R(s);
    seed(R);
  }

  void mwc128x32::seed(result_type x, result_type y, result_type z, result_type c) {
    S.x = x;
    S.y = y;
    S.z = z;
    S.c = clamp_carry(c);
  }

  // Equality comparable concept
  bool operator==(const mwc128x32 &R1, const mwc128x32 &R2) {
    return R1.S == R2.S;
  }

  bool operator!=(const mwc128x32 &R1, const mwc128x32 &R2) {
    return not(R1 == R2);
  }

  // Parallel random number generator concept

  // Other useful methods
  const char *const mwc128x32::name_str = "mwc128x32";

  const char *mwc128x32::name() {
    return name_str;
  }

}  // namespace trng
