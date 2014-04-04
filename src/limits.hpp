// Copyright (c) 2000-2014, Heiko Bauke
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
    
    //using ::std::numeric_limits;
    template<typename T> 
    class numeric_limits {
    public:
      static const bool is_specialized=::std::numeric_limits<T>::is_specialized;
      static T min() throw() {
	return ::std::numeric_limits<T>::min();
      }
      static T max() throw() {
	return ::std::numeric_limits<T>::max();
      }
      static const int  digits=::std::numeric_limits<T>::digits;
      static const int  digits10=::std::numeric_limits<T>::digits10;
      static const bool is_signed=::std::numeric_limits<T>::is_signed;
      static const bool is_integer=::std::numeric_limits<T>::is_integer;
      static const bool is_exact=::std::numeric_limits<T>::is_exact;
      static const int radix=::std::numeric_limits<T>::radix;
      static T epsilon() throw() {
	return ::std::numeric_limits<T>::epsilon();
      }
      static T round_error() throw() {
	return ::std::numeric_limits<T>::round_error();
      }      
      static const int  min_exponent=::std::numeric_limits<T>::min_exponent;
      static const int  min_exponent10=::std::numeric_limits<T>::min_exponent10;
      static const int  max_exponent=::std::numeric_limits<T>::max_exponent;
      static const int  max_exponent10=::std::numeric_limits<T>::max_exponent10;
      static const bool has_infinity=::std::numeric_limits<T>::has_infinity;
      static const bool has_quiet_NaN=::std::numeric_limits<T>::has_quiet_NaN;
      static const bool has_signaling_NaN=::std::numeric_limits<T>::has_signaling_NaN;
      static const ::std::float_denorm_style has_denorm=::std::numeric_limits<T>::has_denorm;
      static const bool has_denorm_loss=::std::numeric_limits<T>::has_denorm_loss;
      static T infinity() throw() {
	return ::std::numeric_limits<T>::infinity();
      }
      static T quiet_NaN() throw() {
	return ::std::numeric_limits<T>::quiet_NaN();
      }
      static T signaling_NaN() throw() {
	return ::std::numeric_limits<T>::signaling_NaN();
      }
      static T denorm_min() throw() {
	return ::std::numeric_limits<T>::nenom_min();
      }
      static const bool is_iec559=::std::numeric_limits<T>::is_iec559;
      static const bool is_bounded=::std::numeric_limits<T>::is_bounded;
      static const bool is_modulo=::std::numeric_limits<T>::is_modulo;
      static const bool traps=::std::numeric_limits<T>::traps;
      static const bool tinyness_before=::std::numeric_limits<T>::tinyness_before;
      static const ::std::float_round_style round_style=::std::numeric_limits<T>::round_style;
    };

    // --- float ---

    template<> 
    class numeric_limits<float> {
    public:
      static const bool is_specialized=::std::numeric_limits<float>::is_specialized;
      TRNG_CUDA_ENABLE
      static float min() throw() {
#if defined __CUDA_ARCH__
	return FLT_MIN;
#else
	return ::std::numeric_limits<float>::min();
#endif
      }
      TRNG_CUDA_ENABLE
      static float max() throw() {
#if defined __CUDA_ARCH__
	return CUDART_MAX_NORMAL_F;
#else
	return ::std::numeric_limits<float>::max();
#endif
      }
      static const int  digits=::std::numeric_limits<float>::digits;
      static const int  digits10=::std::numeric_limits<float>::digits10;
      static const bool is_signed=::std::numeric_limits<float>::is_signed;
      static const bool is_integer=::std::numeric_limits<float>::is_integer;
      static const bool is_exact=::std::numeric_limits<float>::is_exact;
      static const int radix=::std::numeric_limits<float>::radix;
      TRNG_CUDA_ENABLE
      static float epsilon() throw() {
#if defined __CUDA_ARCH__
	return FLT_EPSILON;
#else
	return ::std::numeric_limits<float>::epsilon();
#endif
      }
      TRNG_CUDA_ENABLE
      static float round_error() throw() {
#if defined __CUDA_ARCH__
	return 0.5f;
#else
	return ::std::numeric_limits<float>::round_error();
#endif
      }      
      static const int  min_exponent=::std::numeric_limits<float>::min_exponent;
      static const int  min_exponent10=::std::numeric_limits<float>::min_exponent10;
      static const int  max_exponent=::std::numeric_limits<float>::max_exponent;
      static const int  max_exponent10=::std::numeric_limits<float>::max_exponent10;
      static const bool has_infinity=::std::numeric_limits<float>::has_infinity;
      static const bool has_quiet_NaN=::std::numeric_limits<float>::has_quiet_NaN;
      static const bool has_signaling_NaN=::std::numeric_limits<float>::has_signaling_NaN;
      static const ::std::float_denorm_style has_denorm=::std::numeric_limits<float>::has_denorm;
      static const bool has_denorm_loss=::std::numeric_limits<float>::has_denorm_loss;
      TRNG_CUDA_ENABLE
      static float infinity() throw() {
#if defined __CUDA_ARCH__
	return CUDART_INF_F;
#else
	return ::std::numeric_limits<float>::infinity();
#endif
      }
      TRNG_CUDA_ENABLE
      static float quiet_NaN() throw() {
#if defined __CUDA_ARCH__
	return CUDART_NAN_F;
#else
	return ::std::numeric_limits<float>::quiet_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static float signaling_NaN() throw() {
#if defined __CUDA_ARCH__
	return CUDART_NAN_F;
#else
	return ::std::numeric_limits<float>::signaling_NaN();
#endif
      }
      TRNG_CUDA_ENABLE
      static float denorm_min() throw() {
#if defined __CUDA_ARCH__
	return CUDART_MIN_DENORM_F;
#else
	return ::std::numeric_limits<float>::denorm_min();
#endif
      }
      static const bool is_iec559=::std::numeric_limits<float>::is_iec559;
      static const bool is_bounded=::std::numeric_limits<float>::is_bounded;
      static const bool is_modulo=::std::numeric_limits<float>::is_modulo;
      static const bool traps=::std::numeric_limits<float>::traps;
      static const bool tinyness_before=::std::numeric_limits<float>::tinyness_before;
      static const ::std::float_round_style round_style=::std::numeric_limits<float>::round_style;
    };
    
    // --- double ---

    template<> 
    class numeric_limits<double> {
    public:
      static const bool is_specialized=::std::numeric_limits<double>::is_specialized;
      TRNG_CUDA_ENABLE
      static double min() throw() {
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
      static double max() throw() {
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
      static const int  digits=::std::numeric_limits<double>::digits;
      static const int  digits10=::std::numeric_limits<double>::digits10;
      static const bool is_signed=::std::numeric_limits<double>::is_signed;
      static const bool is_integer=::std::numeric_limits<double>::is_integer;
      static const bool is_exact=::std::numeric_limits<double>::is_exact;
      static const int radix=::std::numeric_limits<double>::radix;
      TRNG_CUDA_ENABLE
      static double epsilon() throw() {
#if defined __CUDA_ARCH__
	return DBL_EPSILON;
#else
	return ::std::numeric_limits<double>::epsilon();
#endif
      }
      TRNG_CUDA_ENABLE
      static double round_error() throw() {
#if defined __CUDA_ARCH__
	return 0.5;
#else
	return ::std::numeric_limits<double>::round_error();
#endif
      }      
      static const int  min_exponent=::std::numeric_limits<double>::min_exponent;
      static const int  min_exponent10=::std::numeric_limits<double>::min_exponent10;
      static const int  max_exponent=::std::numeric_limits<double>::max_exponent;
      static const int  max_exponent10=::std::numeric_limits<double>::max_exponent10;
      static const bool has_infinity=::std::numeric_limits<double>::has_infinity;
      static const bool has_quiet_NaN=::std::numeric_limits<double>::has_quiet_NaN;
      static const bool has_signaling_NaN=::std::numeric_limits<double>::has_signaling_NaN;
      static const ::std::float_denorm_style has_denorm=::std::numeric_limits<double>::has_denorm;
      static const bool has_denorm_loss=::std::numeric_limits<double>::has_denorm_loss;
      TRNG_CUDA_ENABLE
      static double infinity() throw() {
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
      static double quiet_NaN() throw() {
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
      static double signaling_NaN() throw() {
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
      static double denorm_min() throw() {
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
      static const bool is_iec559=::std::numeric_limits<double>::is_iec559;
      static const bool is_bounded=::std::numeric_limits<double>::is_bounded;
      static const bool is_modulo=::std::numeric_limits<double>::is_modulo;
      static const bool traps=::std::numeric_limits<double>::traps;
      static const bool tinyness_before=::std::numeric_limits<double>::tinyness_before;
      static const ::std::float_round_style round_style=::std::numeric_limits<double>::round_style;
    };
    
  }

}

#endif
