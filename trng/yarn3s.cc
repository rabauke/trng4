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

#include "yarn3s.hpp"

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes
  const yarn3s::parameter_type yarn3s::trng0 =
      parameter_type(2025213985, 1112953677, 2038969601);
  const yarn3s::parameter_type yarn3s::trng1 = parameter_type(1287767370, 1045931779, 58150106);

  // Random number engine concept
  yarn3s::yarn3s(yarn3s::parameter_type P) : P{P} {}

  yarn3s::yarn3s(unsigned long s, yarn3s::parameter_type P) : P{P} { seed(s); }

  void yarn3s::seed() { (*this) = yarn3s(); }

  void yarn3s::seed(unsigned long s) {
    int64_t t(s);
    t %= modulus;
    if (t < 0)
      t += modulus;
    S.r[0] = static_cast<result_type>(t);
    S.r[1] = 1;
    S.r[2] = 1;
  }

  void yarn3s::seed(yarn3s::result_type s1, yarn3s::result_type s2, yarn3s::result_type s3) {
    S.r[0] = s1 % modulus;
    if (S.r[0] < 0)
      S.r[0] += modulus;
    S.r[1] = s2 % modulus;
    if (S.r[1] < 0)
      S.r[1] += modulus;
    S.r[2] = s3 % modulus;
    if (S.r[2] < 0)
      S.r[2] += modulus;
  }

  // Equality comparable concept
  bool operator==(const yarn3s &R1, const yarn3s &R2) { return R1.P == R2.P and R1.S == R2.S; }

  bool operator!=(const yarn3s &R1, const yarn3s &R2) { return not(R1 == R2); }

  // other useful methods
  const char *const yarn3s::name_str = "yarn3s";

  const char *yarn3s::name() { return name_str; }

  const int_math::power<yarn3s::modulus, yarn3s::gen> yarn3s::g;

}  // namespace trng
