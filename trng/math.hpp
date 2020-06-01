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

#if !(defined TRNG_MATH_HPP)

#define TRNG_MATH_HPP

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
#include <cmath>
#else
#include <mathimf.h>
#include <cmath>
#endif

#include <trng/cuda.hpp>
#include <ciso646>

namespace trng {

  namespace math {

    // --- trigonometric functions ---------------------------------------

    using ::std::sin;
    using ::std::cos;
    using ::std::tan;

    inline float sec(float x) { return 1.0f / cos(x); }
    inline double sec(double x) { return 1.0 / cos(x); }
    inline long double sec(long double x) { return 1.0l / cos(x); }

    inline float csc(float x) { return 1.0f / sin(x); }
    inline double csc(double x) { return 1.0 / sin(x); }
    inline long double csc(long double x) { return 1.0l / sin(x); }

#if defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL)
    inline float cot(float x) { return ::cotf(x); }
    using ::cot;
    inline long double cot(long double x) { return ::cotl(x); }
#else
    inline float cot(float x) { return cos(x) / sin(x); }
    inline double cot(double x) { return cos(x) / sin(x); }
    inline long double cot(long double x) { return cos(x) / sin(x); }
#endif

    using ::std::asin;
    using ::std::acos;
    using ::std::atan;
    using ::std::atan2;

    inline float atan(float y, float x) { return atan2(y, x); }
    inline double atan(double y, double x) { return atan2(y, x); }
    inline long double atan(long double y, long double x) { return atan2(y, x); }

    // --- hyperbolic functions ------------------------------------------

    using ::std::sinh;
    using ::std::cosh;
    using ::std::tanh;

    inline float sech(float x) { return 1.0f / cosh(x); }
    inline double sech(double x) { return 1.0 / cosh(x); }
    inline long double sech(long double x) { return 1.0l / cosh(x); }

    inline float csch(float x) { return 1.0f / sinh(x); }
    inline double csch(double x) { return 1.0 / sinh(x); }
    inline long double csch(long double x) { return 1.0l / sinh(x); }

    inline float coth(float x) { return cosh(x) / sinh(x); }
    inline double coth(double x) { return cosh(x) / sinh(x); }
    inline long double coth(long double x) { return cosh(x) / sinh(x); }

    using ::std::asinh;
    using ::std::acosh;
    using ::std::atanh;

    inline float asech(float x) { return ::std::acosh(1.0f / x); }
    inline double asech(double x) { return ::std::acosh(1.0 / x); }
    inline long double asech(long double x) { return ::std::acosh(1.0l / x); }

    inline float acsch(float x) { return ::std::asinh(1.0f / x); }
    inline double acsch(double x) { return ::std::asinh(1.0 / x); }
    inline long double acsch(long double x) { return ::std::asinh(1.0l / x); }

    inline float acoth(float x) { return 0.5f * ::std::log((x + 1.0f) / (x - 1.0f)); }
    inline double acoth(double x) { return 0.5 * ::std::log((x + 1.0) / (x - 1.0)); }
    inline long double acoth(long double x) {
      return 0.5l * ::std::log((x + 1.0l) / (x - 1.0l));
    }

    // --- exponential functions -----------------------------------------

    using ::std::exp;

#if defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL)
    inline float exp10(float x) { return ::exp10f(x); }
    using ::exp10;
    inline long double exp10(long double x) { return ::exp10l(x); }
#else
    inline float exp10(float x) { return ::std::pow(10.0f, x); }
    inline double exp10(double x) { return ::std::pow(10.0, x); }
    inline long double exp10(long double x) { return ::std::pow(10.0l, x); }
#endif

    using ::std::exp2;
    using ::std::expm1;

    using ::std::frexp;
    using ::std::ldexp;
    using ::std::log;

    TRNG_CUDA_ENABLE
    inline float ln(float x) { return log(x); }

    TRNG_CUDA_ENABLE
    inline double ln(double x) { return log(x); }

    inline long double ln(long double x) { return log(x); }

    inline float log(float b, float x) { return log(x) * log(b); }
    inline double log(double b, double x) { return log(x) * log(b); }
    inline long double log(long double b, long double x) { return log(x) * log(b); }

    using ::std::log10;
    using ::std::log2;

    TRNG_CUDA_ENABLE
    inline float ln1p(float x) { return ::std::log1p(x); }

    TRNG_CUDA_ENABLE
    inline double ln1p(double x) { return ::std::log1p(x); }

    inline long double ln1p(long double x) { return ::std::log1p(x); }

    using ::std::pow;
    using ::std::sqrt;
    inline float sqrt(float r, float x) { return ::std::pow(x, 1.0f / r); }
    inline double sqrt(double r, double x) { return ::std::pow(x, 1.0 / r); }
    inline long double sqrt(long double r, long double x) { return ::std::pow(x, 1.0l / r); }

    // --- nearest integer functions -------------------------------------

    using ::std::ceil;
    using ::std::floor;
    using ::std::modf;
    using ::std::round;

    // --- misc functions ------------------------------------------------

    using ::std::fmod;
    using ::std::abs;
    using ::std::fabs;

    inline float frac(float x) {
      float dummy;
      return ::std::modf(x, &dummy);
    }
    inline double frac(double x) {
      double dummy;
      return ::std::modf(x, &dummy);
    }
    inline long double frac(long double x) {
      long double dummy;
      return ::std::modf(x, &dummy);
    }

    using ::std::isfinite;
    using ::std::isnan;
    using ::std::isnormal;

  }  // namespace math

}  // namespace trng

#endif
