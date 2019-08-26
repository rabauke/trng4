// Copyright (c) 2000-2019, Heiko Bauke
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

#if !(defined TRNG_SPECIAL_FUNCTIONS_HPP)

#define TRNG_SPECIAL_FUNCTIONS_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/math.hpp>
#include <trng/constants.hpp>
#include <trng/utility.hpp>
#include <cerrno>
#include <algorithm>
#include <ciso646>

namespace trng {

  namespace math {

    TRNG_CUDA_ENABLE
    inline float ln_Gamma(float x) { return ::std::lgamma(x); }

    TRNG_CUDA_ENABLE
    inline double ln_Gamma(double x) { return ::std::lgamma(x); }

#if !(defined __CUDA_ARCH__)
    inline long double ln_Gamma(long double x) { return ::std::lgamma(x); }
#endif

    TRNG_CUDA_ENABLE
    inline float Gamma(float x) { return ::std::tgamma(x); }

    TRNG_CUDA_ENABLE
    inline double Gamma(double x) { return ::std::tgamma(x); }

#if !(defined __CUDA_ARCH__)
    inline long double Gamma(long double x) { return ::std::tgamma(x); }
#endif

    // --- Beta function -----------------------------------------------

    TRNG_CUDA_ENABLE
    inline float Beta(float x, float y) {
      return exp(ln_Gamma(x) + ln_Gamma(y) - ln_Gamma(x + y));
    }

    TRNG_CUDA_ENABLE
    inline double Beta(double x, double y) {
      return exp(ln_Gamma(x) + ln_Gamma(y) - ln_Gamma(x + y));
    }

#if !(defined __CUDA_ARCH__)
    inline long double Beta(long double x, long double y) {
      return exp(ln_Gamma(x) + ln_Gamma(y) - ln_Gamma(x + y));
    }
#endif

    // --- ln of binomial coefficient ----------------------------------

    TRNG_CUDA_ENABLE
    inline float ln_binomial(float n, float m) {
      return ln_Gamma(n + 1.f) - ln_Gamma(m + 1.f) - ln_Gamma(n - m + 1.f);
    }

    TRNG_CUDA_ENABLE
    inline double ln_binomial(double n, double m) {
      return ln_Gamma(n + 1.) - ln_Gamma(m + 1.) - ln_Gamma(n - m + 1.);
    }

#if !(defined __CUDA_ARCH__)
    inline long double ln_binomial(long double n, long double m) {
      return ln_Gamma(n + 1.l) - ln_Gamma(m + 1.l) - ln_Gamma(n - m + 1.l);
    }
#endif

    // --- Pochhammer function -----------------------------------------

    inline float Pochhammer(float x, float a) { return exp(ln_Gamma(x + a) - ln_Gamma(x)); }

    inline double Pochhammer(double x, double a) { return exp(ln_Gamma(x + a) - ln_Gamma(x)); }

#if !(defined __CUDA_ARCH__)
    inline long double Pochhammer(long double x, long double a) {
      return exp(ln_Gamma(x + a) - ln_Gamma(x));
    }
#endif

    // --- incomplete Gamma functions ----------------------------------

    namespace detail {

      // compute incomplete Gamma function
      //
      //  gamma(a, x) = Int(exp(-t)*t^(a-1), t=0..x)
      //
      // or
      //
      //  P(a, x) = gamma(a, x) / Gamma(a)
      //
      // by series expansion
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaP_ser(T a, T x) {
        const int itmax = 32;
        const T eps = T(4) * numeric_limits<T>::epsilon();
        if (x < eps)
          return T(0);
        T xx(T(1) / a), n(a), sum(xx);
        int i(0);
        do {
          ++n;
          ++i;
          xx *= x / n;
          sum += xx;
        } while (abs(xx) > eps * abs(sum) and i < itmax);
        if (by_Gamma_a)
          return exp(-x + a * ln(x) - math::ln_Gamma(a)) * sum;
        return exp(-x + a * ln(x)) * sum;
      }

      // compute complementary incomplete Gamma function
      //
      //  Gamma(a, x) = Int(exp(-t)*t^(a-1), t=x..oo)
      //
      // or
      //
      //  Q(a, x) = Gamma(a, x) / Gamma(a) = 1 - P(a, x)
      //
      // by continued fraction
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaQ_cf(T a, T x) {
        const T itmax = T(32);
        const T eps = T(4) * numeric_limits<T>::epsilon();
        const T min = T(4) * numeric_limits<T>::min();
        // set up for evaluating continued fraction by modied Lentz's method
        T del, bi(x + T(1) - a), ci(T(1) / min), di(T(1) / bi), h(di), i(T(0));
        do {  // iterate
          ++i;
          T ai = -i * (i - a);
          bi += T(2);
          di = ai * di + bi;
          if (abs(di) < min)
            di = min;
          ci = bi + ai / ci;
          if (abs(ci) < min)
            ci = min;
          di = T(1) / di;
          del = di * ci;
          h *= del;
        } while ((abs(del - T(1)) > eps) and i < itmax);
        if (by_Gamma_a)
          return exp(-x + a * ln(x) - math::ln_Gamma(a)) * h;
        return exp(-x + a * ln(x)) * h;
      }

      // P(a, x) and gamma(a, x)
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaP(T a, T x) {
        if (x < T(0) or a <= T(0))
          return numeric_limits<T>::signaling_NaN();
        if (by_Gamma_a) {
          if (x < a + T(1))
            return GammaP_ser<T, true>(a, x);
          return T(1) - GammaQ_cf<T, true>(a, x);
        }
        if (x < a + T(1))
          return GammaP_ser<T, false>(a, x);
        return math::Gamma(a) - GammaQ_cf<T, false>(a, x);
      }

      // Q(a, x) and Gamma(a, x)
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaQ(T a, T x) {
        if (x < T(0) or a <= T(0))
          return numeric_limits<T>::signaling_NaN();
        if (by_Gamma_a) {
          if (x < a + T(1))
            return T(1) - GammaP_ser<T, true>(a, x);
          return GammaQ_cf<T, true>(a, x);
        }
        if (x < a + T(1))
          return math::Gamma(a) - GammaP_ser<T, false>(a, x);
        return GammaQ_cf<T, false>(a, x);
      }

    }  // namespace detail

    // P(x, a)
    TRNG_CUDA_ENABLE
    inline float GammaP(float a, float x) { return detail::GammaP<float, true>(a, x); }

    TRNG_CUDA_ENABLE
    inline double GammaP(double a, double x) { return detail::GammaP<double, true>(a, x); }

#if !(defined __CUDA_ARCH__)
    inline long double GammaP(long double a, long double x) {
      return detail::GammaP<long double, true>(a, x);
    }
#endif

    // Q(x, a)
    TRNG_CUDA_ENABLE
    inline float GammaQ(float a, float x) { return detail::GammaQ<float, true>(a, x); }

    TRNG_CUDA_ENABLE
    inline double GammaQ(double a, double x) { return detail::GammaQ<double, true>(a, x); }

#if !(defined __CUDA_ARCH__)
    inline long double GammaQ(long double a, long double x) {
      return detail::GammaQ<long double, true>(a, x);
    }
#endif

    // gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float inc_gamma(float a, float x) { return detail::GammaP<float, false>(a, x); }

    TRNG_CUDA_ENABLE
    inline double inc_gamma(double a, double x) { return detail::GammaP<double, false>(a, x); }

#if !(defined __CUDA_ARCH__)
    inline long double inc_gamma(long double a, long double x) {
      return detail::GammaP<long double, false>(a, x);
    }
#endif

    // Gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float cinc_gamma(float a, float x) { return detail::GammaQ<float, false>(a, x); }

    TRNG_CUDA_ENABLE
    inline double cinc_gamma(double a, double x) { return detail::GammaQ<double, false>(a, x); }

#if !(defined __CUDA_ARCH__)
    inline long double cinc_gamma(long double a, long double x) {
      return detail::GammaQ<long double, false>(a, x);
    }
#endif

    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE T inv_GammaP(T a, T p) {
        const T eps = sqrt(numeric_limits<T>::epsilon()), a1 = a - T(1),
                glna = math::ln_Gamma(a), lna1 = ln(a1), afac = exp(a1 * (lna1 - T(1)) - glna);
        T x, t;
        // initial guess
        if (a > T(1)) {
          const T pp = p < T(1) / T(2) ? p : T(1) - p;
          t = sqrt(-T(2) * ln(pp));
          x = static_cast<T>((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t);
          x = p < T(1) / T(2) ? -x : x;
          x = utility::max(T(1) / T(1000),
                           a * pow(T(1) - T(1) / (T(9) * a) - x / (T(3) * sqrt(a)), T(3)));
        } else {
          t = static_cast<T>(1.0 - a * (0.253 + a * 0.12));
          x = p < t ? (pow(p / t, T(1) / a)) : (T(1) - ln(T(1) - (p - t) / (T(1) - t)));
        }
        // refinement by Halley's method
        for (int i = 0; i < 16; ++i) {
          if (x < T(0)) {
            x = T(0);
            break;
          }
          const T err = GammaP<T, true>(a, x) - p;
          if (a > T(1))
            t = afac * exp(-(x - a1) + a1 * (ln(x) - lna1));
          else
            t = exp(-x + a1 * ln(x) - glna);
          const T u = err / t;
          t = u / (T(1) - utility::min(T(1), u * ((a - T(1)) / x - T(1))) / T(2));
          x -= t;
          x = x <= T(0) ? (x + t) / T(2) : x;
          if (abs(t) < eps * x)
            break;
        }
        return x;
      }

    }  // namespace detail

    // inverse of GammaP
    TRNG_CUDA_ENABLE
    inline float inv_GammaP(float a, float p) { return detail::inv_GammaP(a, p); }

    // inverse of GammaP
    TRNG_CUDA_ENABLE
    inline double inv_GammaP(double a, double p) { return detail::inv_GammaP(a, p); }

    // inverse of GammaP
#if !(defined __CUDA_ARCH__)
    inline long double inv_GammaP(long double a, long double p) {
      return detail::inv_GammaP(a, p);
    }
#endif

    // --- regularized incomplete Beta function ------------------------

    // see Applied Statistics (1973), vol.22, no.3, pp.409--411
    // algorithm AS 63
    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE T Beta_I(T x, T p, T q, T norm) {
        if (p <= 0 or q <= 0 or x < 0 or x > 1) {
#if !(defined __CUDA_ARCH__)
          errno = EDOM;
#endif
          return numeric_limits<T>::quiet_NaN();
        }
        const T eps = 4 * numeric_limits<T>::epsilon();
        T psq = p + q, cx = 1 - x;
        const bool flag = (p < psq * x);
        if (flag) {
          // use  I(x, p, q) = 1-I(1-x, q, p)
          utility::swap(x, cx);
          utility::swap(p, q);
        }
        T term = 1, i = 1, y = 1, rx = x / cx, temp = q - i;
        int s = static_cast<int>(q + cx * psq);
        if (s == 0)
          rx = x;
        while (true) {
          term *= temp * rx / (p + i);
          y += term;
          temp = abs(term);
          if (temp <= eps and temp <= eps * y)
            break;
          i++;
          s--;
          if (s >= 0) {
            temp = q - i;
            if (s == 0)
              rx = x;
          } else {
            temp = psq;
            psq++;
          }
        }
        y *= exp(p * ln(x) + (q - 1) * ln(cx)) / p / norm;
        return flag ? 1 - y : y;
      }

    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float Beta_I(float x, float p, float q, float norm) {
      return detail::Beta_I(x, p, q, norm);
    }

    TRNG_CUDA_ENABLE
    inline float Beta_I(float x, float p, float q) {
      return detail::Beta_I(x, p, q, Beta(p, q));
    }

    TRNG_CUDA_ENABLE
    inline double Beta_I(double x, double p, double q, double norm) {
      return detail::Beta_I(x, p, q, norm);
    }

    TRNG_CUDA_ENABLE
    inline double Beta_I(double x, double p, double q) {
      return detail::Beta_I(x, p, q, Beta(p, q));
    }

#if !(defined __CUDA_ARCH__)
    inline long double Beta_I(long double x, long double p, long double q, long double norm) {
      return detail::Beta_I(x, p, q, norm);
    }

    inline long double Beta_I(long double x, long double p, long double q) {
      return detail::Beta_I(x, p, q, Beta(p, q));
    }
#endif

    // --- inverse of regularized incomplete Beta function -------------

    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE inline T inv_Beta_I(T x, T p, T q, T norm) {
        if (x < math::numeric_limits<T>::epsilon())
          return 0;
        if (1 - x < math::numeric_limits<T>::epsilon())
          return 1;
        // solve via Newton method
        T y(T(1) / T(2));
        if (2 * p >= 1 and 2 * q >= 1)
          y = (3 * p - 1) / (3 * p + 3 * q - 2);  // the approximate median
        for (int i = 0; i < math::numeric_limits<T>::digits; ++i) {
          const T f(math::Beta_I(y, p, q, norm) - x);
          const T df(math::pow(1 - y, q - 1) * math::pow(y, p - 1) / norm);
          T dy(f / df);
          if (math::abs(dy) < 2 * math::numeric_limits<T>::epsilon())
            break;
          // avoid overshooting
          while (y - dy < 0 or y - dy > 1) {
            dy *= 3;
            dy /= 4;
          }
          y -= dy;
        }
        return y;
      }

    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float inv_Beta_I(float x, float p, float q, float norm) {
      return detail::inv_Beta_I(x, p, q, norm);
    }

    TRNG_CUDA_ENABLE
    inline float inv_Beta_I(float x, float p, float q) {
      return detail::inv_Beta_I(x, p, q, Beta(p, q));
    }

    TRNG_CUDA_ENABLE
    inline double inv_Beta_I(double x, double p, double q, double norm) {
      return detail::inv_Beta_I(x, p, q, norm);
    }

    TRNG_CUDA_ENABLE
    inline double inv_Beta_I(double x, double p, double q) {
      return detail::inv_Beta_I(x, p, q, Beta(p, q));
    }

#if !(defined __CUDA_ARCH__)
    inline long double inv_Beta_I(long double x, long double p, long double q,
                                  long double norm) {
      return detail::inv_Beta_I(x, p, q, norm);
    }

    inline long double inv_Beta_I(long double x, long double p, long double q) {
      return detail::inv_Beta_I(x, p, q, Beta(p, q));
    }
#endif

    // --- error function ----------------------------------------------

    TRNG_CUDA_ENABLE
    inline float erf(float x) { return ::std::erf(x); }

    TRNG_CUDA_ENABLE
    inline double erf(double x) { return ::std::erf(x); }

#if !(defined __CUDA_ARCH__)
    inline long double erf(long double x) { return ::std::erf(x); }
#endif

    // --- complementary error function --------------------------------

    TRNG_CUDA_ENABLE
    inline float erfc(float x) { return ::std::erfc(x); }

    TRNG_CUDA_ENABLE
    inline double erfc(double x) { return ::std::erfc(x); }

#if !(defined __CUDA_ARCH__)
    inline long double erfc(long double x) { return std::erfc(x); }
#endif

    // --- normal distribution function  -------------------------------

    TRNG_CUDA_ENABLE
    inline float Phi(float x) {
      return 0.5f + 0.5f * erf(constants<float>::one_over_sqrt_2() * x);
    }

    TRNG_CUDA_ENABLE
    inline double Phi(double x) {
      return 0.5 + 0.5 * erf(constants<double>::one_over_sqrt_2() * x);
    }

    inline long double Phi(long double x) {
      return 0.5l + 0.5l * erf(constants<long double>::one_over_sqrt_2() * x);
    }

    // --- inverse of normal distribution function  --------------------

    // this function is based on an approximation by Peter J. Acklam
    // see http://home.online.no/~pjacklam/notes/invnorm/ for details

    namespace detail {

      template<typename T>
      struct inv_Phi_traits;

      template<>
      struct inv_Phi_traits<float> {
        TRNG_CUDA_ENABLE
        static float a(int i) noexcept {
          const float a_[] = {-3.969683028665376e+01f, 2.209460984245205e+02f,
                              -2.759285104469687e+02f, 1.383577518672690e+02f,
                              -3.066479806614716e+01f, 2.506628277459239e+00f};
          return a_[i];
        }
        TRNG_CUDA_ENABLE
        static float b(int i) noexcept {
          const float b_[] = {-5.447609879822406e+01f, 1.615858368580409e+02f,
                              -1.556989798598866e+02f, 6.680131188771972e+01f,
                              -1.328068155288572e+01f};
          return b_[i];
        }
        TRNG_CUDA_ENABLE
        static float c(int i) noexcept {
          const float c_[] = {-7.784894002430293e-03f, -3.223964580411365e-01f,
                              -2.400758277161838e+00f, -2.549732539343734e+00f,
                              4.374664141464968e+00f,  2.938163982698783e+00f};
          return c_[i];
        }
        TRNG_CUDA_ENABLE
        static float d(int i) noexcept {
          const float d_[] = {7.784695709041462e-03f, 3.224671290700398e-01f,
                              2.445134137142996e+00f, 3.754408661907416e+00f};
          return d_[i];
        }
        TRNG_CUDA_ENABLE
        static float x_low() noexcept { return 0.02425f; }
        TRNG_CUDA_ENABLE
        static float x_high() noexcept { return 1.0f - 0.02425f; }
        TRNG_CUDA_ENABLE
        static float zero() noexcept { return 0.0f; }
        TRNG_CUDA_ENABLE
        static float one() noexcept { return 1.0f; }
        TRNG_CUDA_ENABLE
        static float one_half() noexcept { return 0.5f; }
        TRNG_CUDA_ENABLE
        static float minus_two() noexcept { return -2.0f; }
      };

      template<>
      struct inv_Phi_traits<double> {
        TRNG_CUDA_ENABLE
        static double a(int i) noexcept {
          const double a_[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                               -2.759285104469687e+02, 1.383577518672690e+02,
                               -3.066479806614716e+01, 2.506628277459239e+00};
          return a_[i];
        }
        TRNG_CUDA_ENABLE
        static double b(int i) noexcept {
          const double b_[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                               -1.556989798598866e+02, 6.680131188771972e+01,
                               -1.328068155288572e+01};
          return b_[i];
        }
        TRNG_CUDA_ENABLE
        static double c(int i) noexcept {
          const double c_[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                               -2.400758277161838e+00, -2.549732539343734e+00,
                               4.374664141464968e+00,  2.938163982698783e+00};
          return c_[i];
        }
        TRNG_CUDA_ENABLE
        static double d(int i) noexcept {
          const double d_[] = {7.784695709041462e-03, 3.224671290700398e-01,
                               2.445134137142996e+00, 3.754408661907416e+00};
          return d_[i];
        }
        TRNG_CUDA_ENABLE
        static double x_low() noexcept { return 0.02425; }
        TRNG_CUDA_ENABLE
        static double x_high() noexcept { return 1.0 - 0.02425; }
        TRNG_CUDA_ENABLE
        static double zero() noexcept { return 0.0; }
        TRNG_CUDA_ENABLE
        static double one() noexcept { return 1.0; }
        TRNG_CUDA_ENABLE
        static double one_half() noexcept { return 0.5; }
        TRNG_CUDA_ENABLE
        static double minus_two() noexcept { return -2.0; }
      };

      template<>
      struct inv_Phi_traits<long double> {
        TRNG_CUDA_ENABLE
        static long double a(int i) noexcept {
          const long double a_[] = {-3.969683028665376e+01l, 2.209460984245205e+02l,
                                    -2.759285104469687e+02l, 1.383577518672690e+02l,
                                    -3.066479806614716e+01l, 2.506628277459239e+00l};
          return a_[i];
        }
        TRNG_CUDA_ENABLE
        static long double b(int i) noexcept {
          const long double b_[] = {-5.447609879822406e+01l, 1.615858368580409e+02l,
                                    -1.556989798598866e+02l, 6.680131188771972e+01l,
                                    -1.328068155288572e+01l};
          return b_[i];
        }
        TRNG_CUDA_ENABLE
        static long double c(int i) noexcept {
          const long double c_[] = {-7.784894002430293e-03l, -3.223964580411365e-01l,
                                    -2.400758277161838e+00l, -2.549732539343734e+00l,
                                    4.374664141464968e+00l,  2.938163982698783e+00l};
          return c_[i];
        }
        TRNG_CUDA_ENABLE
        static long double d(int i) noexcept {
          const long double d_[] = {7.784695709041462e-03l, 3.224671290700398e-01l,
                                    2.445134137142996e+00l, 3.754408661907416e+00l};
          return d_[i];
        }
        TRNG_CUDA_ENABLE
        static long double x_low() noexcept { return 0.02425l; }
        TRNG_CUDA_ENABLE
        static long double x_high() noexcept { return 1.0l - 0.02425l; }
        TRNG_CUDA_ENABLE
        static long double zero() noexcept { return 0.0l; }
        TRNG_CUDA_ENABLE
        static long double one() noexcept { return 1.0l; }
        TRNG_CUDA_ENABLE
        static long double one_half() noexcept { return 0.5l; }
        TRNG_CUDA_ENABLE
        static long double minus_two() noexcept { return -2.0l; }
      };

      template<typename T>
      TRNG_CUDA_ENABLE T inv_Phi(T x) {
        if (x < inv_Phi_traits<T>::zero() or x > inv_Phi_traits<T>::one()) {
#if !(defined __CUDA_ARCH__)
          errno = EDOM;
#endif
          return numeric_limits<T>::quiet_NaN();
        }
        if (x == inv_Phi_traits<T>::zero())
          return -numeric_limits<T>::infinity();
        if (x == inv_Phi_traits<T>::one())
          return numeric_limits<T>::infinity();
        T t, q;
        if (x < inv_Phi_traits<T>::x_low()) {
          // Rational approximation for lower region
          q = sqrt(inv_Phi_traits<T>::minus_two() * ln(x));
          t = (((((inv_Phi_traits<T>::c(0) * q + inv_Phi_traits<T>::c(1)) * q +
                  inv_Phi_traits<T>::c(2)) *
                     q +
                 inv_Phi_traits<T>::c(3)) *
                    q +
                inv_Phi_traits<T>::c(4)) *
                   q +
               inv_Phi_traits<T>::c(5)) /
              ((((inv_Phi_traits<T>::d(0) * q + inv_Phi_traits<T>::d(1)) * q +
                 inv_Phi_traits<T>::d(2)) *
                    q +
                inv_Phi_traits<T>::d(3)) *
                   q +
               inv_Phi_traits<T>::one());
        } else if (x < inv_Phi_traits<T>::x_high()) {
          // Rational approximation for central region
          q = x - inv_Phi_traits<T>::one_half();
          T r = q * q;
          t = (((((inv_Phi_traits<T>::a(0) * r + inv_Phi_traits<T>::a(1)) * r +
                  inv_Phi_traits<T>::a(2)) *
                     r +
                 inv_Phi_traits<T>::a(3)) *
                    r +
                inv_Phi_traits<T>::a(4)) *
                   r +
               inv_Phi_traits<T>::a(5)) *
              q /
              (((((inv_Phi_traits<T>::b(0) * r + inv_Phi_traits<T>::b(1)) * r +
                  inv_Phi_traits<T>::b(2)) *
                     r +
                 inv_Phi_traits<T>::b(3)) *
                    r +
                inv_Phi_traits<T>::b(4)) *
                   r +
               inv_Phi_traits<T>::one());
        } else {
          // Rational approximation for upper region
          q = sqrt(inv_Phi_traits<T>::minus_two() * ln(1 - x));
          t = -(((((inv_Phi_traits<T>::c(0) * q + inv_Phi_traits<T>::c(1)) * q +
                   inv_Phi_traits<T>::c(2)) *
                      q +
                  inv_Phi_traits<T>::c(3)) *
                     q +
                 inv_Phi_traits<T>::c(4)) *
                    q +
                inv_Phi_traits<T>::c(5)) /
              ((((inv_Phi_traits<T>::d(0) * q + inv_Phi_traits<T>::d(1)) * q +
                 inv_Phi_traits<T>::d(2)) *
                    q +
                inv_Phi_traits<T>::d(3)) *
                   q +
               inv_Phi_traits<T>::one());
        }
        // refinement by Halley rational method
        if (numeric_limits<T>::epsilon() < 1e-9) {
          T e(Phi(t) - x);
          T u(e * constants<T>::sqrt_2pi() * exp(t * t * inv_Phi_traits<T>::one_half()));
          t -= u / (inv_Phi_traits<T>::one() + t * u * inv_Phi_traits<T>::one_half());
        }
        return t;
      }

    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float inv_Phi(float x) { return detail::inv_Phi<float>(x); }

    TRNG_CUDA_ENABLE
    inline double inv_Phi(double x) { return detail::inv_Phi<double>(x); }

#if !(defined __CUDA_ARCH__)
    inline long double inv_Phi(long double x) { return detail::inv_Phi<long double>(x); }
#endif

    // --- inverse of error function  ----------------------------------

    // see http://mathworld.wolfram.com/InverseErf.html
    // see The On-Line Encyclopedia of Integer Sequences!
    // http://www.research.att.com/~njas/sequences/A007019
    // http://www.research.att.com/~njas/sequences/A092676

    TRNG_CUDA_ENABLE
    inline float inv_erf(float x) {
      if (abs(x) < 1.0f / 8.0f) {
        x *= 0.886226925452758013649085f;  // sqrt(pi)/2
        float x2 = x * x, x3 = x2 * x, x4 = x2 * x2;
        return x + (1.0f / 3.0f + 7.0f / 30.0f * x2 + 127.0f / 630.0f * x4) * x3;
      }
      return inv_Phi(0.5f * (x + 1.0f)) * constants<float>::one_over_sqrt_2();
    }

    TRNG_CUDA_ENABLE
    inline double inv_erf(double x) {
      if (abs(x) < 1.0 / 20.0) {
        x *= 0.886226925452758013649085;  // sqrt(pi)/2
        double x2 = x * x, x3 = x2 * x, x4 = x2 * x2, x5 = x3 * x2;
        return x + (1.0 / 3.0 + 127.0 / 630.0 * x4) * x3 +
               (7.0 / 30.0 + 4369.0 / 22680.0 * x4) * x5;
      }
      return inv_Phi(0.5 * (x + 1.0)) * constants<double>::one_over_sqrt_2();
    }

#if !(defined __CUDA_ARCH__)
    inline long double inv_erf(long double x) {
      if (abs(x) < 1.0l / 24.0l) {
        x *= 0.886226925452758013649085l;  // sqrt(pi)/2
        long double x2 = x * x, x3 = x2 * x, x4 = x2 * x2, x5 = x3 * x2, x7 = x3 * x4;
        return x + 1.0l / 3.0l * x3 + 7.0l / 30.0l * x5 + 127.0l / 630.0l * x7 +
               4369.0l / 22680.0l * x4 * x5 + 34807.0l / 178200.0l * x4 * x7;
      }
      return inv_Phi(0.5l * (x + 1.0l)) * constants<long double>::one_over_sqrt_2();
    }
#endif

    // --- inverse of complementary error function  --------------------

    TRNG_CUDA_ENABLE
    inline float inv_erfc(float x) {
      return -inv_Phi(0.5f * x) * constants<float>::one_over_sqrt_2();
    }

    TRNG_CUDA_ENABLE
    inline double inv_erfc(double x) {
      return -inv_Phi(0.5 * x) * constants<double>::one_over_sqrt_2();
    }

#if !(defined __CUDA_ARCH__)
    inline long double inv_erfc(long double x) {
      return -inv_Phi(0.5l * x) * constants<long double>::one_over_sqrt_2();
    }
#endif

  }  // namespace math

}  // namespace trng

#endif
