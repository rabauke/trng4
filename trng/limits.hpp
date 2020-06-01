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

#if !(defined TRNG_LIMITS_HPP)

#define TRNG_LIMITS_HPP

#include <limits>
#include <cfloat>
#include <trng/cuda.hpp>

#if defined TRNG_CUDA
#include <math_constants.h>
#endif

namespace trng {

  namespace math {

    // using ::std::numeric_limits;
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

    // --- float ---

    template<>
    class numeric_limits<float> {
    public:
      static constexpr bool is_specialized = ::std::numeric_limits<float>::is_specialized;
      TRNG_CUDA_ENABLE
      static constexpr float min() noexcept {
#if defined __CUDA_ARCH__
        return FLT_MIN;
#else
        return ::std::numeric_limits<float>::min();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr float max() noexcept {
#if defined __CUDA_ARCH__
        return CUDART_MAX_NORMAL_F;
#else
        return ::std::numeric_limits<float>::max();
#endif
      }
      static constexpr int digits = ::std::numeric_limits<float>::digits;
      static constexpr int digits10 = ::std::numeric_limits<float>::digits10;
      static constexpr bool is_signed = ::std::numeric_limits<float>::is_signed;
      static constexpr bool is_integer = ::std::numeric_limits<float>::is_integer;
      static constexpr bool is_exact = ::std::numeric_limits<float>::is_exact;
      static constexpr int radix = ::std::numeric_limits<float>::radix;
      TRNG_CUDA_ENABLE
      static constexpr float epsilon() noexcept {
#if defined __CUDA_ARCH__
        return FLT_EPSILON;
#else
        return ::std::numeric_limits<float>::epsilon();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr float round_error() noexcept {
#if defined __CUDA_ARCH__
        return 0.5f;
#else
        return ::std::numeric_limits<float>::round_error();
#endif
      }
      static constexpr int min_exponent = ::std::numeric_limits<float>::min_exponent;
      static constexpr int min_exponent10 = ::std::numeric_limits<float>::min_exponent10;
      static constexpr int max_exponent = ::std::numeric_limits<float>::max_exponent;
      static constexpr int max_exponent10 = ::std::numeric_limits<float>::max_exponent10;
      static constexpr bool has_infinity = ::std::numeric_limits<float>::has_infinity;
      static constexpr bool has_quiet_NaN = ::std::numeric_limits<float>::has_quiet_NaN;
      static constexpr bool has_signaling_NaN = ::std::numeric_limits<float>::has_signaling_NaN;
      static constexpr ::std::float_denorm_style has_denorm =
          ::std::numeric_limits<float>::has_denorm;
      static constexpr bool has_denorm_loss = ::std::numeric_limits<float>::has_denorm_loss;
      TRNG_CUDA_ENABLE
      static constexpr float infinity() noexcept {
#if defined __CUDA_ARCH__
        return CUDART_INF_F;
#else
        return ::std::numeric_limits<float>::infinity();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr float quiet_NaN() noexcept {
#if defined __CUDA_ARCH__
        return CUDART_NAN_F;
#else
        return ::std::numeric_limits<float>::quiet_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr float signaling_NaN() noexcept {
#if defined __CUDA_ARCH__
        return CUDART_NAN_F;
#else
        return ::std::numeric_limits<float>::signaling_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr float denorm_min() noexcept {
#if defined __CUDA_ARCH__
        return CUDART_MIN_DENORM_F;
#else
        return ::std::numeric_limits<float>::denorm_min();
#endif
      }
      static constexpr bool is_iec559 = ::std::numeric_limits<float>::is_iec559;
      static constexpr bool is_bounded = ::std::numeric_limits<float>::is_bounded;
      static constexpr bool is_modulo = ::std::numeric_limits<float>::is_modulo;
      static constexpr bool traps = ::std::numeric_limits<float>::traps;
      static constexpr bool tinyness_before = ::std::numeric_limits<float>::tinyness_before;
      static constexpr ::std::float_round_style round_style =
          ::std::numeric_limits<float>::round_style;
    };

    // --- double ---

    template<>
    class numeric_limits<double> {
    public:
      static constexpr bool is_specialized = ::std::numeric_limits<double>::is_specialized;
      TRNG_CUDA_ENABLE
      static constexpr double min() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return DBL_MIN;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::min();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr double max() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return DBL_MAX;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::max();
#endif
      }
      static constexpr int digits = ::std::numeric_limits<double>::digits;
      static constexpr int digits10 = ::std::numeric_limits<double>::digits10;
      static constexpr bool is_signed = ::std::numeric_limits<double>::is_signed;
      static constexpr bool is_integer = ::std::numeric_limits<double>::is_integer;
      static constexpr bool is_exact = ::std::numeric_limits<double>::is_exact;
      static constexpr int radix = ::std::numeric_limits<double>::radix;
      TRNG_CUDA_ENABLE
      static constexpr double epsilon() noexcept {
#if defined __CUDA_ARCH__
        return DBL_EPSILON;
#else
        return ::std::numeric_limits<double>::epsilon();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr double round_error() noexcept {
#if defined __CUDA_ARCH__
        return 0.5;
#else
        return ::std::numeric_limits<double>::round_error();
#endif
      }
      static constexpr int min_exponent = ::std::numeric_limits<double>::min_exponent;
      static constexpr int min_exponent10 = ::std::numeric_limits<double>::min_exponent10;
      static constexpr int max_exponent = ::std::numeric_limits<double>::max_exponent;
      static constexpr int max_exponent10 = ::std::numeric_limits<double>::max_exponent10;
      static constexpr bool has_infinity = ::std::numeric_limits<double>::has_infinity;
      static constexpr bool has_quiet_NaN = ::std::numeric_limits<double>::has_quiet_NaN;
      static constexpr bool has_signaling_NaN =
          ::std::numeric_limits<double>::has_signaling_NaN;
      static constexpr ::std::float_denorm_style has_denorm =
          ::std::numeric_limits<double>::has_denorm;
      static constexpr bool has_denorm_loss = ::std::numeric_limits<double>::has_denorm_loss;
      TRNG_CUDA_ENABLE
      static constexpr double infinity() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return CUDART_INF;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::infinity();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr double quiet_NaN() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return CUDART_NAN;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::quiet_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr double signaling_NaN() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return CUDART_NAN;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::signaling_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static constexpr double denorm_min() noexcept {
#if defined __CUDA_ARCH__
#if __CUDA_ARCH__ >= 130
        return CUDART_MIN_DENORM;
#else
        return 0;  // make compiler happy
#endif
#else
        return ::std::numeric_limits<double>::denorm_min();
#endif
      }
      static constexpr bool is_iec559 = ::std::numeric_limits<double>::is_iec559;
      static constexpr bool is_bounded = ::std::numeric_limits<double>::is_bounded;
      static constexpr bool is_modulo = ::std::numeric_limits<double>::is_modulo;
      static constexpr bool traps = ::std::numeric_limits<double>::traps;
      static constexpr bool tinyness_before = ::std::numeric_limits<double>::tinyness_before;
      static constexpr ::std::float_round_style round_style =
          ::std::numeric_limits<double>::round_style;
    };

  }  // namespace math

}  // namespace trng

#endif
