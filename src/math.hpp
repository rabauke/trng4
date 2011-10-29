#if !(defined TRNG_MATH_HPP)

#define TRNG_MATH_HPP

#if !(defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL))
#  include <cmath>
#else
#  include <mathimf.h>
#  include <cmath>
#endif

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
      return log(x+::std::sqrt(x*x-1.0f));
    }
    inline double acosh(double x) {
      return log(x+::std::sqrt(x*x-1.0));
    }
    inline long double acosh(long double x) {
      return log(x+::std::sqrt(x*x-1.0l));
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
    
    inline float ln(float x) {
      return log(x);
    }

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
      return ::std::pow(x, 1.0f/x);
    }
    inline double sqrt(double r, double x) {
      return ::std::pow(x, 1.0/x);
    }
    inline long double sqrt(long double r, long double x) {
      return ::std::pow(x, 1.0l/x);
    }
  
    // --- nearest integer functions -------------------------------------

    using ::std::ceil;
    using ::std::floor;
  
    inline long lrint(float x) {
      return ::lrintf(x);
    }
    using ::lrint;
    inline long lrint(long double x) {
      return ::lrintl(x);
    }

    inline long long llrint(float x) {
      return ::llrintf(x);
    }
    using ::llrint;
    inline long long llrint(long double x) {
      return ::llrintl(x);
    }
  
    inline long lround(float x) {
      return ::lroundf(x);
    }
    using ::lround;
    inline long lround(long double x) {
      return ::lroundl(x);
    }

    inline long long llround(float x) {
      return ::llroundf(x);
    }
    using ::llround;
    inline long long llround(long double x) {
      return ::llroundl(x);
    }

    using ::std::modf;

    inline float nearbyint(float x) {
      return ::nearbyintf(x);
    }
    using ::nearbyint;
    inline long double nearbyint(long double x) {
      return ::nearbyintl(x);
    }

    inline float rint(float x) {
      return ::rintf(x);
    }
    using ::rint;
    inline long double rint(long double x) {
      return ::rintl(x);
    }

    inline float round(float x) {
      return ::roundf(x);
    }
    using ::round;
    inline long double round(long double x) {
      return ::roundl(x);
    }

    inline float trunc(float x) {
      return ::truncf(x);
    }
    using ::trunc;
    inline long double trunc(long double x) {
      return ::truncl(x);
    }

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
