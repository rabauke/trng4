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

#if !(defined TRNG_LIMITS_HPP)

#define TRNG_LIMITS_HPP

#include <limits>
#include <cfloat>
#include <trng/cuda.hpp>

#if defined __CUDACC__ && !(defined __HIPCC__)
#include <math_constants.h>
#include <cuda/std/limits>
#endif

namespace trng {

  namespace math {

#if defined __CUDACC__ && !(defined __HIPCC__)
    using cuda::std::numeric_limits;
#elif defined __HIPCC__
    template<typename T>
    class numeric_limits {
    public:
      static constexpr bool is_specialized = ::std::numeric_limits<T>::is_specialized;
      static constexpr T min() noexcept { return ::std::numeric_limits<T>::min(); }
      static constexpr T max() noexcept { return ::std::numeric_limits<T>::max(); }
      static constexpr int digits = ::std::numeric_limits<T>::digits;
      static constexpr int digits10 = ::std::numeric_limits<T>::digits10;
      static constexpr bool is_signed = ::std::numeric_limits<T>::is_signed;
      static constexpr bool is_integer = ::std::numeric_limits<T>::is_integer;
      static constexpr bool is_exact = ::std::numeric_limits<T>::is_exact;
      static constexpr int radix = ::std::numeric_limits<T>::radix;
      static constexpr T epsilon() noexcept { return ::std::numeric_limits<T>::epsilon(); }
      static constexpr T round_error() noexcept {
        return ::std::numeric_limits<T>::round_error();
      }
      static constexpr int min_exponent = ::std::numeric_limits<T>::min_exponent;
      static constexpr int min_exponent10 = ::std::numeric_limits<T>::min_exponent10;
      static constexpr int max_exponent = ::std::numeric_limits<T>::max_exponent;
      static constexpr int max_exponent10 = ::std::numeric_limits<T>::max_exponent10;
      static constexpr bool has_infinity = ::std::numeric_limits<T>::has_infinity;
      static constexpr bool has_quiet_NaN = ::std::numeric_limits<T>::has_quiet_NaN;
      static constexpr bool has_signaling_NaN = ::std::numeric_limits<T>::has_signaling_NaN;
      static constexpr ::std::float_denorm_style has_denorm =
          ::std::numeric_limits<T>::has_denorm;
      static constexpr bool has_denorm_loss = ::std::numeric_limits<T>::has_denorm_loss;
      static constexpr T infinity() noexcept { return ::std::numeric_limits<T>::infinity(); }
      static constexpr T quiet_NaN() noexcept { return ::std::numeric_limits<T>::quiet_NaN(); }
      static constexpr T signaling_NaN() noexcept {
        return ::std::numeric_limits<T>::signaling_NaN();
      }
      static constexpr T denorm_min() noexcept {
        return ::std::numeric_limits<T>::denorm_min();
      }
      static constexpr bool is_iec559 = ::std::numeric_limits<T>::is_iec559;
      static constexpr bool is_bounded = ::std::numeric_limits<T>::is_bounded;
      static constexpr bool is_modulo = ::std::numeric_limits<T>::is_modulo;
      static constexpr bool traps = ::std::numeric_limits<T>::traps;
      static constexpr bool tinyness_before = ::std::numeric_limits<T>::tinyness_before;
      static constexpr ::std::float_round_style round_style =
          ::std::numeric_limits<T>::round_style;
    };
#else
    using std::numeric_limits;
#endif

  }  // namespace math

}  // namespace trng

#endif
