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

#if !(defined TRNG_MATH_HPP)

#define TRNG_MATH_HPP

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
#  include <cmath>
#else
#  include <mathimf.h>
#  include <cmath>
#endif

#include <trng/cuda.hpp>

namespace trng {
  
  namespace math {
  
    // --- trigonometric functions ---------------------------------------

    using ::std::sin;
    using ::std::cos;
    using ::std::tan;

    inline float sec(float x) {
      return 1.0f/cos(x);
    }
    inline double sec(double x) {
      return 1.0/cos(x);
    }
    inline long double sec(long double x) {
      return 1.0l/cos(x);
    }

    inline float csc(float x) {
      return 1.0f/sin(x);
    }
    inline double csc(double x) {
      return 1.0/sin(x);
    }
    inline long double csc(long double x) {
      return 1.0l/sin(x);
    }

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float cot(float x) {
      return cos(x)/sin(x);
    }
    inline double cot(double x) {
      return cos(x)/sin(x);
    }
    inline long double cot(long double x) {
      return cos(x)/sin(x);
    }
#else
    inline float cot(float x) {
      return ::cotf(x);
    }
    using ::cot;
    inline long double cot(long double x) {
      return ::cotl(x);
    }
#endif

    using ::std::asin;
    using ::std::acos;
    using ::std::atan;
    using ::std::atan2;

    inline float atan(float y, float x) {
      return atan2(y, x);
    }
    inline double atan(double y, double x) {
      return atan2(y, x);
    }
    inline long double atan(long double y, long double x) {
      return atan2(y, x);
    }

    // --- hyperbolic functions ------------------------------------------

    using ::std::sinh;
    using ::std::cosh;
    using ::std::tanh;

    inline float sech(float x) {
      return 1.0f/cosh(x);
    }
    inline double sech(double x) {
      return 1.0/cosh(x);
    }
    inline long double sech(long double x) {
      return 1.0l/cosh(x);
    }

    inline float csch(float x) {
      return 1.0f/sinh(x);
    }
    inline double csch(double x) {
      return 1.0/sinh(x);
    }
    inline long double csch(long double x) {
      return 1.0l/sinh(x);
    }

    inline float coth(float x) {
      return cosh(x)/sinh(x);
    }
    inline double coth(double x) {
      return cosh(x)/sinh(x);
    }
    inline long double coth(long double x) {
      return cosh(x)/sinh(x);
    }

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float asinh(float x) {
      return ::std::log(x+::std::sqrt(x*x+1.0f));
    }
    inline double asinh(double x) {
      return ::std::log(x+::std::sqrt(x*x+1.0));
    }
    inline long double asinh(long double x) {
      return ::std::log(x+::std::sqrt(x*x+1.0l));
    }
#else
    inline float asinh(float x) {
      return ::asinhf(x);
    }
    using ::asinh;
    inline long double asinh(long double x) {
      return ::asinhl(x);
    }
#endif

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float acosh(float x) {
      return ::std::log(x+::std::sqrt(x*x-1.0f));
    }
    inline double acosh(double x) {
      return ::std::log(x+::std::sqrt(x*x-1.0));
    }
    inline long double acosh(long double x) {
      return ::std::log(x+::std::sqrt(x*x-1.0l));
    }
#else
    inline float acosh(float x) {
      return ::acoshf(x);
    }
    using ::acosh;
    inline long double acosh(long double x) {
      return ::acoshl(x);
    }
#endif

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float atanh(float x) {
      return 0.5f*::std::log((1.0f+x)/(1.0f-x));
    }
    inline double atanh(double x) {
      return 0.5*::std::log((1.0+x)/(1.0-x));
    }
    inline long double atanh(long double x) {
      return 0.5l*::std::log((1.0l+x)/(1.0l-x));
    }
#else
    inline float atanh(float x) {
      return ::atanhf(x);
    }
    using ::atanh;
    inline long double atanh(long double x) {
      return ::atanhl(x);
    }
#endif

    inline float asech(float x) {
      return ::std::log((1.0f+::std::sqrt(1.0f-x*x))/x);
    }
    inline double asech(double x) {
      return ::std::log((1.0+::std::sqrt(1.0-x*x))/x);
    }
    inline long double asech(long double x) {
      return ::std::log((1.0l+::std::sqrt(1.0l-x*x))/x);
    }
  
    inline float acsch(float x) {
      float t=1.0f/x;
      return ::std::log(t+::std::sqrt(1.0f+t*t));
    }
    inline double acsch(double x) {
      double t=1.0/x;
      return ::std::log(t+::std::sqrt(1.0+t*t));
    }
    inline long double acsch(long double x) {
      long double t=1.0l/x;
      return ::std::log(t+::std::sqrt(1.0l+t*t));
    }

    inline float acoth(float x) {
      return 0.5f*::std::log((x+1.0f)/(x-1.0f));
    }
    inline double acoth(double x) {
      return 0.5*::std::log((x+1.0)/(x-1.0));
    }
    inline long double acoth(long double x) {
      return 0.5l*::std::log((x+1.0l)/(x-1.0l));
    }

    // --- exponential functions -----------------------------------------

    using ::std::exp;

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float exp10(float x) {
      return ::std::pow(10.0f, x);
    }
    inline double exp10(double x) {
      return ::std::pow(10.0, x);
    }
    inline long double exp10(long double x) {
      return ::std::pow(10.0l, x);
    }
#else
    inline float exp10(float x) {
      return ::exp10f(x);
    }
    using ::exp10;
    inline long double exp10(long double x) {
      return ::exp10l(x);
    }
#endif

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float exp2(float x) {
      return ::std::pow(2.0f, x);
    }
    inline double exp2(double x) {
      return ::std::pow(2.0, x);
    }
    inline long double exp2(long double x) {
      return ::std::pow(2.0l, x);
    }
#else
    inline float exp2(float x) {
      return ::exp2f(x);
    }
    using ::exp2;
    inline long double exp2(long double x) {
      return ::exp2l(x);
    }
#endif

    using ::std::frexp;
    using ::std::ldexp;
    using ::std::log;

    template<typename T>
    inline T log2_floor(T x) {
      T y(0);
      while (x>0) {
	x>>=1;
	++y;
      };
      --y;
      return y;
    }
    
    template<typename T>
    inline T log2_ceil(T x) {
      T y(log2_floor(x));
      if ((T(1)<<y)<x)
        ++y;
      return y;
    }
    
    template<typename T>
    inline T pow2(T x) {
      return T(1) << x;
    }
    
    TRNG_CUDA_ENABLE
    inline float ln(float x) {
      return log(x);
    }

    TRNG_CUDA_ENABLE
    inline double ln(double x) {
      return log(x);
    }

    inline long double ln(long double x) {
      return log(x);
    }

    inline float log(float b, float x) {
      static float last_b(2.0f);
      static float last_log_b(1.0f/log(last_b));
      if (b!=last_b) {
	last_b=b;
	last_log_b=1.0f/log(b);
      }
      return log(x)*last_log_b;
    }
    inline double log(double b, double x) {
      static double last_b(2.0);
      static double last_log_b(1.0/log(last_b));
      if (b!=last_b) {
	last_b=b;
	last_log_b=1.0/log(b);
      }
      return log(x)*last_log_b;
    }
    inline long double log(long double b, long double x) {
      static long double last_b(2.0l);
      static long double last_log_b(1.0l/log(last_b));
      if (b!=last_b) {
	last_b=b;
	last_log_b=1.0l/log(b);
      }
      return log(x)*last_log_b;
    }

    using ::std::log10;

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float log2(float x) {
      return log(x)*1.44269504088896340735992f;
    }
    inline double log2(double x) {
      return log(x)*1.44269504088896340735992;
    }
    inline long double log2(long double x) {
      return log(x)*1.44269504088896340735992l;
    }
#else
    inline float log2(float x) {
      return ::log2f(x);
    }
    using ::log2;
    inline long double log2(long double x) {
      return ::log2l(x);
    }
#endif

    using ::std::pow;
    using ::std::sqrt;
    inline float sqrt(float r, float x) {
      return ::std::pow(x, 1.0f/r);
    }
    inline double sqrt(double r, double x) {
      return ::std::pow(x, 1.0/r);
    }
    inline long double sqrt(long double r, long double x) {
      return ::std::pow(x, 1.0l/r);
    }
  
    // --- nearest integer functions -------------------------------------

    using ::std::ceil;
    using ::std::floor;
    using ::std::modf;

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
    inline float round(float x) {
      return floor(x+0.5f);
    }
    inline double round(double x) {
      return floor(x+0.5);
    }
    inline long double round(long double x) {
      return floor(x+0.5l);
    }
#else
    inline float round(float x) {
      return ::roundf(x);
    }
    using ::round;
    inline long double round(long double x) {
      return ::roundl(x);
    }
#endif

    // --- misc functions ------------------------------------------------

    using ::std::fmod;
    using ::std::abs;
    using ::std::fabs;

    inline float frac(float x) {
      return x-floor(x);
    }
    inline double frac(double x) {
      return x-floor(x);
    }
    inline long double frac(long double x) {
      return x-floor(x);
    }

  }

}

#endif
