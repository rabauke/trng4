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

#include "yarn5s.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes
  const yarn5s::parameter_type yarn5s::trng0 =
      parameter_type(1053223373, 1530818118, 1612122482, 133497989, 573245311);
  const yarn5s::parameter_type yarn5s::trng1 =
      parameter_type(2068619238, 2138332912, 671754166, 1442240992, 1526958817);

  // Random number engine concept
  yarn5s::yarn5s(const parameter_type &P) : P{P} {}

  yarn5s::yarn5s(unsigned long s, const parameter_type &P) : P{P} { seed(s); }

  void yarn5s::seed() { (*this) = yarn5s(); }

  void yarn5s::seed(unsigned long s) {
    int64_t t = s;
    t %= modulus;
    if (t < 0)
      t += modulus;
    S.r[0] = static_cast<result_type>(t);
    S.r[1] = 1;
    S.r[2] = 1;
    S.r[3] = 1;
    S.r[4] = 1;
  }

  void yarn5s::seed(yarn5s::result_type s1, yarn5s::result_type s2, yarn5s::result_type s3,
                    yarn5s::result_type s4, yarn5s::result_type s5) {
    S.r[0] = s1 % modulus;
    if (S.r[0] < 0)
      S.r[0] += modulus;
    S.r[1] = s2 % modulus;
    if (S.r[1] < 0)
      S.r[1] += modulus;
    S.r[2] = s3 % modulus;
    if (S.r[2] < 0)
      S.r[2] += modulus;
    S.r[3] = s4 % modulus;
    if (S.r[3] < 0)
      S.r[3] += modulus;
    S.r[4] = s5 % modulus;
    if (S.r[4] < 0)
      S.r[4] += modulus;
  }

  // Equality comparable concept
  bool operator==(const yarn5s &R1, const yarn5s &R2) { return R1.P == R2.P and R1.S == R2.S; }

  bool operator!=(const yarn5s &R1, const yarn5s &R2) { return not(R1 == R2); }

  // other useful methods
  const char *const yarn5s::name_str = "yarn5s";

  const char *yarn5s::name() { return name_str; }

  const int_math::power<yarn5s::modulus, yarn5s::gen> yarn5s::g;

}  // namespace trng
