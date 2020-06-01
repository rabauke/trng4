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
      return ln_Gamma(n + 1) - ln_Gamma(m + 1) - ln_Gamma(n - m + 1);
    }

    TRNG_CUDA_ENABLE
    inline double ln_binomial(double n, double m) {
      return ln_Gamma(n + 1) - ln_Gamma(m + 1) - ln_Gamma(n - m + 1);
    }

#if !(defined __CUDA_ARCH__)
    inline long double ln_binomial(long double n, long double m) {
      return ln_Gamma(n + 1) - ln_Gamma(m + 1) - ln_Gamma(n - m + 1);
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
      // by series expansion, see "Numerical Recipes" by W. H. Press et al., 3rd edition
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaP_ser(T a, T x) {
        const int itmax{32};
        const T eps{4 * numeric_limits<T>::epsilon()};
        if (x < eps)
          return T{0};
        T xx{1 / a}, n{a}, sum{xx};
        int i{0};
        do {
          ++n;
          ++i;
          xx *= x / n;
          sum += xx;
        } while (abs(xx) > eps * abs(sum) and i < itmax);
#if __cplusplus >= 201703L
        if constexpr (by_Gamma_a)
#else
        if (by_Gamma_a)
#endif
          return exp(-x + a * ln(x) - ln_Gamma(a)) * sum;
        else
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
      // by continued fraction, see "Numerical Recipes" by W. H. Press et al., 3rd edition
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaQ_cf(T a, T x) {
        const T itmax{32};
        const T eps{4 * numeric_limits<T>::epsilon()};
        const T min{4 * numeric_limits<T>::min()};
        // set up for evaluating continued fraction by modified Lentz's method
        T del, bi{x + 1 - a}, ci{1 / min}, di{1 / bi}, h{di}, i{0};
        do {  // iterate
          ++i;
          T ai{-i * (i - a)};
          bi += 2;
          di = ai * di + bi;
          if (abs(di) < min)
            di = min;
          ci = bi + ai / ci;
          if (abs(ci) < min)
            ci = min;
          di = 1 / di;
          del = di * ci;
          h *= del;
        } while ((abs(del - 1) > eps) and i < itmax);
#if __cplusplus >= 201703L
        if constexpr (by_Gamma_a)
#else
        if (by_Gamma_a)
#endif
          return exp(-x + a * ln(x) - ln_Gamma(a)) * h;
        else
          return exp(-x + a * ln(x)) * h;
      }

      // P(a, x) and gamma(a, x)
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaP(T a, T x) {
        if (x < 0 or a <= 0)
          return numeric_limits<T>::signaling_NaN();
#if __cplusplus >= 201703L
        if constexpr (by_Gamma_a) {
#else
        if (by_Gamma_a) {
#endif
          if (x < a + 1)
            return GammaP_ser<T, true>(a, x);
          return 1 - GammaQ_cf<T, true>(a, x);
        } else {
          if (x < a + 1)
            return GammaP_ser<T, false>(a, x);
          return Gamma(a) - GammaQ_cf<T, false>(a, x);
        }
      }

      // Q(a, x) and Gamma(a, x)
      template<typename T, bool by_Gamma_a>
      TRNG_CUDA_ENABLE T GammaQ(T a, T x) {
        if (x < 0 or a <= 0)
          return numeric_limits<T>::signaling_NaN();
#if __cplusplus >= 201703L
        if constexpr (by_Gamma_a) {
#else
        if (by_Gamma_a) {
#endif
          if (x < a + 1)
            return T{1} - GammaP_ser<T, true>(a, x);
          return GammaQ_cf<T, true>(a, x);
        } else {
          if (x < a + 1)
            return Gamma(a) - GammaP_ser<T, false>(a, x);
          return GammaQ_cf<T, false>(a, x);
        }
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
    inline float inc_Gamma(float a, float x) { return detail::GammaQ<float, false>(a, x); }

    TRNG_CUDA_ENABLE
    inline double inc_Gamma(double a, double x) { return detail::GammaQ<double, false>(a, x); }

#if !(defined __CUDA_ARCH__)
    inline long double inc_Gamma(long double a, long double x) {
      return detail::GammaQ<long double, false>(a, x);
    }
#endif

    namespace detail {

      // compute inverse of the incomplete Gamma function p = P(a, x), see "Numerical Recipes"
      // by W. H. Press et al., 3rd edition
      template<typename T>
      TRNG_CUDA_ENABLE T inv_GammaP(T a, T p) {
        const T eps{sqrt(numeric_limits<T>::epsilon())};
        T a1{a - 1};
        T glna{ln_Gamma(a)};
        T lna1{ln(a1)};
        T afac{exp(a1 * (lna1 - 1) - glna)};
        T x;
        // initial guess
        if (a > T{1}) {
          const T pp{p < T{1} / T{2} ? p : 1 - p};
          const T t = {sqrt(-2 * ln(pp))};
          x = (T{2.30753} + t * T{0.27061}) / (1 + t * (T{0.99229} + t * T{0.04481})) - t;
          x = p < T{1} / T{2} ? -x : x;
          x = utility::max(T{1} / T{1000}, a * pow(1 - 1 / (9 * a) - x / (3 * sqrt(a)), T{3}));
        } else {
          const T t{1 - a * (T{0.253} + a * T{0.12})};
          x = p < t ? pow(p / t, 1 / a) : 1 - ln1p(-(p - t) / (1 - t));
        }
        // refinement by Halley's method
        for (int i{0}; i < 16; ++i) {
          if (x <= 0) {
            x = 0;
            break;
          }
          const T err{GammaP<T, true>(a, x) - p};
          T t;
          if (a > 1)
            t = afac * exp(-(x - a1) + a1 * (ln(x) - lna1));
          else
            t = exp(-x + a1 * ln(x) - glna);
          const T u{err / t};
          t = u / (1 - utility::min(T{1}, u * ((a - 1) / x - 1)) / 2);
          x -= t;
          x = x <= 0 ? (x + t) / 2 : x;
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
        const T eps{4 * numeric_limits<T>::epsilon()};
        T psq{p + q}, cx{1 - x};
        const bool flag{p < psq * x};
        if (flag) {
          // use  I(x, p, q) = 1 - I(1 - x, q, p)
          utility::swap(x, cx);
          utility::swap(p, q);
        }
        T term{1}, i{1}, y{1}, rx{x / cx}, temp{q - i};
        int s{static_cast<int>(q + cx * psq)};
        if (s == 0)
          rx = x;
        while (true) {
          term *= temp * rx / (p + i);
          y += term;
          temp = abs(term);
          if (temp <= eps and temp <= eps * y)
            break;
          ++i;
          --s;
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
        if (x < numeric_limits<T>::epsilon())
          return 0;
        if (1 - x < numeric_limits<T>::epsilon())
          return 1;
        // solve via Newton method
        T y{0};
        if (2 * p >= 1 and 2 * q >= 1)
          y = (3 * p - 1) / (3 * p + 3 * q - 2);  // the approximate median
        else {
          // following initial guess given in "Numerical Recipes" by W. H. Press et al., 3rd
          // edition
          const T lnp{ln(p / (p + q))};
          const T lnq{ln(q / (p + q))};
          const T t{exp(p * lnp) / p};
          const T u{exp(q * lnq) / q};
          const T w{t + u};
          if (x < t / w)
            y = pow(p * w * x, 1 / p);
          else
            y = 1 - pow(q * w * (1 - x), 1 / q);
        }
        for (int i{0}; i < numeric_limits<T>::digits; ++i) {
          const T f{Beta_I(y, p, q, norm) - x};
          const T df{pow(1 - y, q - 1) * pow(y, p - 1) / norm};
          T dy(f / df);
          if (abs(f / y) < 2 * numeric_limits<T>::epsilon())
            break;
          // avoid overshooting
          while (y - dy <= 0 or y - dy >= 1)
            dy *= T{3} / T{4};
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

    // --- error function and complementary error function--------------

    using ::std::erf;
    using ::std::erfc;

    // --- normal distribution function  -------------------------------

    TRNG_CUDA_ENABLE
    inline float Phi(float x) {
      x *= constants<float>::one_over_sqrt_2;
      if (x < -0.6744897501960817f * constants<float>::one_over_sqrt_2)
        return 0.5f * erfc(-x);
      if (x > +0.6744897501960817f * constants<float>::one_over_sqrt_2)
        return 1.0f - 0.5f * erfc(x);
      return 0.5f + 0.5f * erf(x);
    }

    TRNG_CUDA_ENABLE
    inline double Phi(double x) {
      x *= constants<double>::one_over_sqrt_2;
      if (x < -0.6744897501960817 * constants<double>::one_over_sqrt_2)
        return 0.5 * erfc(-x);
      if (x > +0.6744897501960817 * constants<double>::one_over_sqrt_2)
        return 1.0 - 0.5 * erfc(x);
      return 0.5 + 0.5 * erf(x);
    }

    inline long double Phi(long double x) {
      x *= constants<long double>::one_over_sqrt_2;
      if (x < -0.6744897501960817l * constants<long double>::one_over_sqrt_2)
        return 0.5 * erfc(-x);
      if (x > +0.6744897501960817l * constants<long double>::one_over_sqrt_2)
        return 1.0l - 0.5l * erfc(x);
      return 0.5l + 0.5l * erf(x);
    }

    // --- inverse of normal distribution function  --------------------

    // this function is based on an approximation by Peter J. Acklam
    // see http://home.online.no/~pjacklam/notes/invnorm/ or
    // https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
    // for details

    namespace detail {

      template<typename T>
      struct inv_Phi_traits {
        static constexpr T a[6]{
            static_cast<T>(-3.969683028665376e+01l), static_cast<T>(2.209460984245205e+02l),
            static_cast<T>(-2.759285104469687e+02l), static_cast<T>(1.383577518672690e+02l),
            static_cast<T>(-3.066479806614716e+01l), static_cast<T>(2.506628277459239e+00l)};
        static constexpr T b[5]{
            static_cast<T>(-5.447609879822406e+01l), static_cast<T>(1.615858368580409e+02l),
            static_cast<T>(-1.556989798598866e+02l), static_cast<T>(6.680131188771972e+01l),
            static_cast<T>(-1.328068155288572e+01l)};
        static constexpr T c[6]{
            static_cast<T>(-7.784894002430293e-03l), static_cast<T>(-3.223964580411365e-01l),
            static_cast<T>(-2.400758277161838e+00l), static_cast<T>(-2.549732539343734e+00l),
            static_cast<T>(4.374664141464968e+00l),  static_cast<T>(2.938163982698783e+00l)};
        static constexpr T d[4]{
            static_cast<T>(7.784695709041462e-03l), static_cast<T>(3.224671290700398e-01l),
            static_cast<T>(2.445134137142996e+00l), static_cast<T>(3.754408661907416e+00l)};
        static constexpr T x_low = static_cast<T>(0.02425l);
        static constexpr T x_high = static_cast<T>(1.0l - 0.02425l);
        static constexpr T one_half = static_cast<T>(0.5l);
      };

      template<typename T>
      constexpr T inv_Phi_traits<T>::a[6];

      template<typename T>
      constexpr T inv_Phi_traits<T>::b[5];

      template<typename T>
      constexpr T inv_Phi_traits<T>::c[6];

      template<typename T>
      constexpr T inv_Phi_traits<T>::d[4];

      // ---------------------------------------------------------------

      template<typename T>
      TRNG_CUDA_ENABLE T inv_Phi_approx(T x) {
        using traits = inv_Phi_traits<T>;
        if (x < 0 or x > 1) {
#if !(defined __CUDA_ARCH__)
          errno = EDOM;
#endif
          return numeric_limits<T>::quiet_NaN();
        }
        if (x == 0)
          return -numeric_limits<T>::infinity();
        if (x == 1)
          return numeric_limits<T>::infinity();
        T t, q;
        if (x < traits::x_low) {
          // Rational approximation for lower region
          q = sqrt(-2 * ln(x));
          t = (((((traits::c[0] * q + traits::c[1]) * q + traits::c[2]) * q + traits::c[3]) *
                    q +
                traits::c[4]) *
                   q +
               traits::c[5]) /
              ((((traits::d[0] * q + traits::d[1]) * q + traits::d[2]) * q + traits::d[3]) * q +
               1);
        } else if (x < traits::x_high) {
          // Rational approximation for central region
          q = x - traits::one_half;
          T r = q * q;
          t = (((((traits::a[0] * r + traits::a[1]) * r + traits::a[2]) * r + traits::a[3]) *
                    r +
                traits::a[4]) *
                   r +
               traits::a[5]) *
              q /
              (((((traits::b[0] * r + traits::b[1]) * r + traits::b[2]) * r + traits::b[3]) *
                    r +
                traits::b[4]) *
                   r +
               1);
        } else {
          // Rational approximation for upper region
          q = sqrt(-2 * ln1p(-x));
          t = -(((((traits::c[0] * q + traits::c[1]) * q + traits::c[2]) * q + traits::c[3]) *
                     q +
                 traits::c[4]) *
                    q +
                traits::c[5]) /
              ((((traits::d[0] * q + traits::d[1]) * q + traits::d[2]) * q + traits::d[3]) * q +
               1);
        }
        return t;
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_Phi(T x) {
        using traits = inv_Phi_traits<T>;
        T y{inv_Phi_approx(x)};
        if (isfinite(y)) {  // refinement by Halley rational method
          const T e{(Phi(y) - x)};
          const T u{e * constants<T>::sqrt_2pi * exp(y * y * traits::one_half)};
          y -= u / (1 + y * u * traits::one_half);
        }
        return y;
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_erf(T x) {
        T y{inv_Phi_approx((x + 1) / 2) * constants<T>::one_over_sqrt_2};
        if (isfinite(y)) {  // refinement by Halley rational method
          const T e{erf(y) - x};
          const T u{e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
          y -= u / (1 + y * u);
        }
        return y;
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_erfc(T x) {
        // step size in the Halley step is proportiaonal to erfc, use symmetry to increase
        // numerical accuracy
        const bool flag{x > 1};
        if (flag)
          x = -(x - 1) + 1;
        T y{-inv_Phi_approx(x / 2) * constants<T>::one_over_sqrt_2};
        if (isfinite(y)) {  // refinement by Halley rational method
          const T e{erfc(y) - x};
          const T u{-e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
          y -= u / (1 + y * u);
        }
        return flag ? -y : y;
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

    TRNG_CUDA_ENABLE
    inline float inv_erf(float x) { return detail::inv_erf(x); }

    TRNG_CUDA_ENABLE
    inline double inv_erf(double x) { return detail::inv_erf(x); }

#if !(defined __CUDA_ARCH__)
    inline long double inv_erf(long double x) { return detail::inv_erf(x); }
#endif

    // --- inverse of complementary error function  --------------------

    TRNG_CUDA_ENABLE
    inline float inv_erfc(float x) { return detail::inv_erfc(x); }

    TRNG_CUDA_ENABLE
    inline double inv_erfc(double x) { return detail::inv_erfc(x); }

#if !(defined __CUDA_ARCH__)
    inline long double inv_erfc(long double x) { return detail::inv_erfc(x); }
#endif

  }  // namespace math

}  // namespace trng

#endif
