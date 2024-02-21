// Copyright (c) 2000-2022, Heiko Bauke
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

    // --- error function and complementary error function--------------

    using ::std::erf;
    using ::std::erfc;

    // --- x - ln(1 + x) -----------------------------------------------

    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE T mln1p(T x) {
        if (abs(x) >= T(1) / T(32))
          return x - ln1p(x);
        // use Taylor expansion for small arguments
        T y{0};
        T x_to_the_n{x * x};
        T sign{1};
        for (int n{2}; n < numeric_limits<T>::digits; ++n) {
          const T delta{sign * x_to_the_n / n};
          y += delta;
          if (abs(delta) < 4 * numeric_limits<T>::epsilon() * y)
            break;
          x_to_the_n *= x;
          sign = -sign;
        }
        return y;
      }

    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float mln1p(float x) {
      return detail::mln1p(x);
    }

    TRNG_CUDA_ENABLE
    inline double mln1p(double x) {
      return detail::mln1p(x);
    }

#if !(defined TRNG_CUDA)
    inline long double mln1p(long double x) {
      return detail::mln1p(x);
    }
#endif

    // --- logarithm of the Gamma function -----------------------------

    TRNG_CUDA_ENABLE
    inline float ln_Gamma(float x) {
      return ::std::lgamma(x);
    }

    TRNG_CUDA_ENABLE
    inline double ln_Gamma(double x) {
      return ::std::lgamma(x);
    }

#if !(defined TRNG_CUDA)
    inline long double ln_Gamma(long double x) {
      return ::std::lgamma(x);
    }
#endif

    // --- Gamma function ----------------------------------------------

    TRNG_CUDA_ENABLE
    inline float Gamma(float x) {
      return ::std::tgamma(x);
    }

    TRNG_CUDA_ENABLE
    inline double Gamma(double x) {
      return ::std::tgamma(x);
    }

#if !(defined TRNG_CUDA)
    inline long double Gamma(long double x) {
      return ::std::tgamma(x);
    }
#endif

    // --- Beta function -----------------------------------------------

    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE T Beta(T x, T y) {
        static const T ln_max{ln(numeric_limits<T>::max())};
        if (x <= 0 or y <= 0) {
#if !(defined TRNG_CUDA)
          errno = EDOM;
#endif
          return numeric_limits<T>::signaling_NaN();
        }
        const T z{x + y};
        if (z * ln(z) - z > ln_max)
          // less accurate but avoids overflow
          return exp(ln_Gamma(x) + ln_Gamma(y) - ln_Gamma(z));
        return Gamma(x) / Gamma(z) * Gamma(y);
      }
    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float Beta(float x, float y) {
      return detail::Beta(x, y);
    }

    TRNG_CUDA_ENABLE
    inline double Beta(double x, double y) {
      return detail::Beta(x, y);
    }

#if !(defined TRNG_CUDA)
    inline long double Beta(long double x, long double y) {
      return detail::Beta(x, y);
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

#if !(defined TRNG_CUDA)
    inline long double ln_binomial(long double n, long double m) {
      return ln_Gamma(n + 1) - ln_Gamma(m + 1) - ln_Gamma(n - m + 1);
    }
#endif

    // --- Pochhammer function -----------------------------------------

    TRNG_CUDA_ENABLE
    inline float Pochhammer(float x, float a) {
      return exp(ln_Gamma(x + a) - ln_Gamma(x));
    }

    TRNG_CUDA_ENABLE
    inline double Pochhammer(double x, double a) {
      return exp(ln_Gamma(x + a) - ln_Gamma(x));
    }

#if !(defined TRNG_CUDA)
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
        const int itmax{numeric_limits<T>::digits};
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
        const T itmax{numeric_limits<T>::digits};
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


      template<typename T>
      class GammaPQ_asympt_coefficients;

      template<>
      class GammaPQ_asympt_coefficients<float> {
      public:
        static constexpr std::size_t n = 15;

        TRNG_CUDA_ENABLE
        float operator[](std::size_t index) const {
          // clang-format off
          static constexpr float d[n]{
            1.f,                                          // 1
            -3.333333333333333333333333333333333333e-01f, // -1 / 3
            8.333333333333333333333333333333333333e-02f,  // 1 / 12
            -1.481481481481481481481481481481481481e-02f, // -2 / 135
            1.157407407407407407407407407407407407e-03f,  // 1 / 864
            3.527336860670194003527336860670194004e-04f,  // 1 / 2835
            -1.787551440329218106995884773662551440e-04f, // -139 / 777600
            3.919263178522437781697040956300215559e-05f,  // 1 / 25515
            -2.185448510679992161473642955124436606e-06f, // -571 / 261273600
            -1.854062210715159960701798836229563253e-06f, // -281 / 151559100
            8.296711340953086005016242131664432272e-07f,  // 163879 / 197522841600
            -1.766595273682607930436005424574240304e-07f, // -5221 / 29554024500
            6.707853543401498580369397100296135722e-09f,  // 5246819 / 782190452736000
            1.026180978424030804257395732272529509e-08f,  // 5459 / 531972441000
            -4.382036018453353186552974622447191234e-09f  // -534703531 / 122021710626816000
          };
          // clang-format on
          return d[index];
        }
      };

      template<>
      class GammaPQ_asympt_coefficients<double> {
      public:
        static constexpr std::size_t n = 27;

        TRNG_CUDA_ENABLE
        double operator[](std::size_t index) const {
          // clang-format off
          static constexpr double d[n]{
            1.,                                          // 1
            -3.333333333333333333333333333333333333e-01, // -1 / 3
            8.333333333333333333333333333333333333e-02,  // 1 / 12
            -1.481481481481481481481481481481481481e-02, // -2 / 135
            1.157407407407407407407407407407407407e-03,  // 1 / 864
            3.527336860670194003527336860670194004e-04,  // 1 / 2835
            -1.787551440329218106995884773662551440e-04, // -139 / 777600
            3.919263178522437781697040956300215559e-05,  // 1 / 25515
            -2.185448510679992161473642955124436606e-06, // -571 / 261273600
            -1.854062210715159960701798836229563253e-06, // -281 / 151559100
            8.296711340953086005016242131664432272e-07,  // 163879 / 197522841600
            -1.766595273682607930436005424574240304e-07, // -5221 / 29554024500
            6.707853543401498580369397100296135722e-09,  // 5246819 / 782190452736000
            1.026180978424030804257395732272529509e-08,  // 5459 / 531972441000
            -4.382036018453353186552974622447191234e-09, // -534703531 / 122021710626816000
            9.147699582236790234182488176331136808e-10,  // 91207079 / 99704934754425000
            -2.551419399494624976687795379938870131e-11, // -4483131259 / 175711263302615040000
            -5.830772132550425067464089450400357975e-11, // -2650986803 / 45465450248017800000
            2.436194802066741624369406967077899429e-11,  // 432261921612371 / 17743323368298066739200000
            -5.027669280114175589090549859257443655e-12, // -6171801683 / 1227567156696480600000
            1.100439203195613477083741744972934113e-13,  // 6232523202521089 / 56636688191607429031526400000
            3.371763262400985378827698841692001848e-13,  // 4283933145517 / 12705320071808574210000000
            -1.392388722418162065919366184895799799e-13, // -25834629665134204969 / 185541790515705937507280486400000
            2.853489380704744320396690990528282989e-14,  // 11963983648109 / 419275562369682948930000000
            -5.139111834242572618990645803004942055e-16, // -1579029138854919086429 / 3072572050940090325120564854784000000
            -1.975228829434944283539624015807109122e-15, // -208697624924077 / 105657441717160103130360000000
            8.099521156704561334071156687025752553e-16   // 746590869962651602203151 / 921771615282027097536169456435200000000
          };
          // clang-format on
          return d[index];
        }
      };

      template<>
      class GammaPQ_asympt_coefficients<long double> {
      public:
        static constexpr std::size_t n = 42;

        long double operator[](std::size_t index) const {
          // clang-format off
          static constexpr long double d[n]{
            1.l,                                          // 1
            -3.333333333333333333333333333333333333e-01l, // -1 / 3
            8.333333333333333333333333333333333333e-02l,  // 1 / 12
            -1.481481481481481481481481481481481481e-02l, // -2 / 135
            1.157407407407407407407407407407407407e-03l,  // 1 / 864
            3.527336860670194003527336860670194004e-04l,  // 1 / 2835
            -1.787551440329218106995884773662551440e-04l, // -139 / 777600
            3.919263178522437781697040956300215559e-05l,  // 1 / 25515
            -2.185448510679992161473642955124436606e-06l, // -571 / 261273600
            -1.854062210715159960701798836229563253e-06l, // -281 / 151559100
            8.296711340953086005016242131664432272e-07l,  // 163879 / 197522841600
            -1.766595273682607930436005424574240304e-07l, // -5221 / 29554024500
            6.707853543401498580369397100296135722e-09l,  // 5246819 / 782190452736000
            1.026180978424030804257395732272529509e-08l,  // 5459 / 531972441000
            -4.382036018453353186552974622447191234e-09l, // -534703531 / 122021710626816000
            9.147699582236790234182488176331136808e-10l,  // 91207079 / 99704934754425000
            -2.551419399494624976687795379938870131e-11l, // -4483131259 / 175711263302615040000
            -5.830772132550425067464089450400357975e-11l, // -2650986803 / 45465450248017800000
            2.436194802066741624369406967077899429e-11l,  // 432261921612371 / 17743323368298066739200000
            -5.027669280114175589090549859257443655e-12l, // -6171801683 / 1227567156696480600000
            1.100439203195613477083741744972934113e-13l,  // 6232523202521089 / 56636688191607429031526400000
            3.371763262400985378827698841692001848e-13l,  // 4283933145517 / 12705320071808574210000000
            -1.392388722418162065919366184895799799e-13l, // -25834629665134204969 / 185541790515705937507280486400000
            2.853489380704744320396690990528282989e-14l,  // 11963983648109 / 419275562369682948930000000
            -5.139111834242572618990645803004942055e-16l, // -1579029138854919086429 / 3072572050940090325120564854784000000
            -1.975228829434944283539624015807109122e-15l, // -208697624924077 / 105657441717160103130360000000
            8.099521156704561334071156687025752553e-16l,  // 746590869962651602203151 / 921771615282027097536169456435200000000
            -1.652253121639816181915148202653511616e-16l, // -29320119130515566117 / 177455371374430493811049182600000000
            2.530543009747888423270610900602673849e-18l,  // 1511513601028097903631961 / 597308006702753559203437807770009600000000
            1.168693973855957658882308765077934746e-17l,  // 2700231121460756431181 / 231046893529508502941986035745200000000
            -4.770037049820484758221678040848165974e-18l, // -8849272268392873147705987190261 / 1855178938018082279529957487152872816640000000000
            9.699126059056237124207096858985853544e-19l,  // 10084288256532215186381 / 10397110208827882632389371608534000000000
            -1.293256553803817501044325677449629120e-20l, // -6208770108287283939483943525987 / 480088045177548944685317694066691261071360000000000
            -6.969230253185693380530546315855194520e-20l, // -6782242429223267933535073 / 97316951554628981439164518255878240000000000
            2.835145432176936599923196187131687508e-20l,  // 2355444393109967510921431436000087153 / 83080196394065199975683597593629056110920990720000000000
            -5.750982159007047500162666139520842346e-21, // -51748587106835353426330148693 / 8998217291595659510809468851493269705120000000000
            6.792953783488914564614965664370431421e-23l,  // 2346608607351903737647919577082115121863 / 34544745660652310149889239879430961530920947941376000000000000
            4.182125426111335857807972609348226866e-22l,  // 7007277101869903281324331583 / 16755301163660883227024528206228847037120000000000
            -1.697153962004760373219505695058247317e-22l, // -2603072187220373277150999431416562396331667 / 15337867073329625706550822506467346919728900885970944000000000000
            3.436215938394319882960425642469428583e-23l,  // 585302872633292617248814587726421 / 17033355386471835584177000251811214753701006400000000000
            -3.643995779628021011969396867845457286e-25l, // -73239727426811935976967471475430268695630993 / 200987410128911415258641978124748114036127517209763250176000000000000
            -2.522535663578433775881227291922899288e-24l  // -110855495796575034381969281033555329 / 43946056897097335807176660649672934064548596512000000000000
          };
          // clang-format on
          return d[index];
        }
      };


      template<typename T>
      TRNG_CUDA_ENABLE T GammaPQ_asympt_R(T a, T eta, T eta_squard_half) {
        GammaPQ_asympt_coefficients<T> coeffs;
        constexpr std::size_t n{coeffs.n - 1};
        T beta[n];
        beta[n - 1] = coeffs[n];
        beta[n - 2] = coeffs[n - 1];
        for (std::ptrdiff_t i{n - 3}; i >= 0; --i)
          beta[i] = beta[i + 2] * (i + 2) / a + coeffs[i + 1];
        T eta_to_the_i{1};
        T y{0};
        for (std::size_t i{0}; i < n; ++i) {
          const T y_old{y};
          y += beta[i] * eta_to_the_i;
          if (y == y_old)
            break;
          eta_to_the_i *= eta;
        }
        return y / (1 + beta[1] / a) * exp(-a * eta_squard_half) / sqrt(a) *
               constants<T>::one_over_sqrt_2pi;
      }


      // calculate regularized incomplete the incomplete Gamma function by an asymptotic
      // expansion, see
      // SIAM Journal on Scientific Computing 34(6) (2012), A2965-A2981
      // https://doi.org/10.48550/arXiv.1306.1754 and references therein
      template<typename T>
      TRNG_CUDA_ENABLE T GammaP_asympt(T a, T x) {
        const T mu{(x - a) / a};
        const T eta_squared_half{mln1p(mu)};
        const T eta{copysign(sqrt(2 * eta_squared_half), mu)};
        const T leading{erfc(-eta * sqrt(a / 2)) / 2};
        const T correction{-GammaPQ_asympt_R(a, eta, eta_squared_half)};
        return leading + correction;
      }


      // calculate regularized complementary incomplete the incomplete Gamma function by an
      // asymptotic expansion, see
      // SIAM Journal on Scientific Computing 34(6) (2012), A2965-A2981
      // https://doi.org/10.48550/arXiv.1306.1754 and references therein
      template<typename T>
      TRNG_CUDA_ENABLE T GammaQ_asympt(T a, T x) {
        const T mu{(x - a) / a};
        const T eta_squared_half{mln1p(mu)};
        const T eta{copysign(sqrt(2 * eta_squared_half), mu)};
        const T leading{erfc(eta * sqrt(a / 2)) / 2};
        const T correction{GammaPQ_asympt_R(a, eta, eta_squared_half)};
        return leading + correction;
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
          if (a > 12 and x > T(3) / T(10) * a and x < T(235) / T(100) * a)
            return GammaP_asympt(a, x);
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
          if (a > 12 and x > T(3) / T(10) * a and x < T(235) / T(100) * a)
            return GammaQ_asympt(a, x);
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
    inline float GammaP(float a, float x) {
      return detail::GammaP<float, true>(a, x);
    }

    TRNG_CUDA_ENABLE
    inline double GammaP(double a, double x) {
      return detail::GammaP<double, true>(a, x);
    }

#if !(defined TRNG_CUDA)
    inline long double GammaP(long double a, long double x) {
      return detail::GammaP<long double, true>(a, x);
    }
#endif

    // Q(x, a)
    TRNG_CUDA_ENABLE
    inline float GammaQ(float a, float x) {
      return detail::GammaQ<float, true>(a, x);
    }

    TRNG_CUDA_ENABLE
    inline double GammaQ(double a, double x) {
      return detail::GammaQ<double, true>(a, x);
    }

#if !(defined TRNG_CUDA)
    inline long double GammaQ(long double a, long double x) {
      return detail::GammaQ<long double, true>(a, x);
    }
#endif

    // gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float inc_gamma(float a, float x) {
      return detail::GammaP<float, false>(a, x);
    }

    TRNG_CUDA_ENABLE
    inline double inc_gamma(double a, double x) {
      return detail::GammaP<double, false>(a, x);
    }

#if !(defined TRNG_CUDA)
    inline long double inc_gamma(long double a, long double x) {
      return detail::GammaP<long double, false>(a, x);
    }
#endif

    // Gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float inc_Gamma(float a, float x) {
      return detail::GammaQ<float, false>(a, x);
    }

    TRNG_CUDA_ENABLE
    inline double inc_Gamma(double a, double x) {
      return detail::GammaQ<double, false>(a, x);
    }

#if !(defined TRNG_CUDA)
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
          const T t{sqrt(-2 * ln(pp))};
          x = (T{2.30753} + t * T{0.27061}) / (1 + t * (T{0.99229} + t * T{0.04481})) - t;
          x = p < T{1} / T{2} ? -x : x;
          x = utility::max(T{1} / T{1000}, a * pow(1 - 1 / (9 * a) - x / (3 * sqrt(a)), T{3}));
        } else {
          const T t{1 - a * (T{0.253} + a * T{0.12})};
          x = p < t ? pow(p / t, 1 / a) : 1 - ln1p(-(p - t) / (1 - t));
        }
        // refinement by Halley's method
        for (int i{0}; i < 32; ++i) {
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
    inline float inv_GammaP(float a, float p) {
      return detail::inv_GammaP(a, p);
    }

    // inverse of GammaP
    TRNG_CUDA_ENABLE
    inline double inv_GammaP(double a, double p) {
      return detail::inv_GammaP(a, p);
    }

    // inverse of GammaP
#if !(defined TRNG_CUDA)
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
#if !(defined TRNG_CUDA)
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

#if !(defined TRNG_CUDA)
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
      TRNG_CUDA_ENABLE T inv_Beta_I(T x, T p, T q, T norm) {
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

#if !(defined TRNG_CUDA)
    inline long double inv_Beta_I(long double x, long double p, long double q,
                                  long double norm) {
      return detail::inv_Beta_I(x, p, q, norm);
    }

    inline long double inv_Beta_I(long double x, long double p, long double q) {
      return detail::inv_Beta_I(x, p, q, Beta(p, q));
    }
#endif

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

#if !(defined TRNG_CUDA)
    inline long double Phi(long double x) {
      x *= constants<long double>::one_over_sqrt_2;
      if (x < -0.6744897501960817l * constants<long double>::one_over_sqrt_2)
        return 0.5l * erfc(-x);
      if (x > +0.6744897501960817l * constants<long double>::one_over_sqrt_2)
        return 1.0l - 0.5l * erfc(x);
      return 0.5l + 0.5l * erf(x);
    }
#endif

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
#if !(defined TRNG_CUDA)
          errno = EDOM;
#endif
          return numeric_limits<T>::quiet_NaN();
        }
        if (x == 0)
          return -numeric_limits<T>::infinity();
        if (x == 1)
          return numeric_limits<T>::infinity();
        if (x < traits::x_low) {
          // rational approximation for lower region
          const T q{sqrt(-2 * ln(x))};
          return (((((traits::c[0] * q + traits::c[1]) * q + traits::c[2]) * q + traits::c[3]) *
                       q +
                   traits::c[4]) *
                      q +
                  traits::c[5]) /
                 ((((traits::d[0] * q + traits::d[1]) * q + traits::d[2]) * q + traits::d[3]) *
                      q +
                  1);
        } else if (x < traits::x_high) {
          // rational approximation for central region
          const T q{x - traits::one_half};
          const T r{q * q};
          return (((((traits::a[0] * r + traits::a[1]) * r + traits::a[2]) * r + traits::a[3]) *
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
          // rational approximation for upper region
          const T q{sqrt(-2 * ln1p(-x))};
          return -(((((traits::c[0] * q + traits::c[1]) * q + traits::c[2]) * q +
                     traits::c[3]) *
                        q +
                    traits::c[4]) *
                       q +
                   traits::c[5]) /
                 ((((traits::d[0] * q + traits::d[1]) * q + traits::d[2]) * q + traits::d[3]) *
                      q +
                  1);
        }
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_Phi(T x) {
        using traits = inv_Phi_traits<T>;
        T y{inv_Phi_approx(x)};
        // refinement by Halley rational method
        if (isfinite(y)) {
          const T e{(Phi(y) - x)};
          const T u{e * constants<T>::sqrt_2pi * exp(y * y * traits::one_half)};
          y -= u / (1 + y * u * traits::one_half);
        }
        // T is a floating point number type with more than 80 bits, a 2nd iteration is
        // required to reach full machine precision
#if __cplusplus >= 201703L
        if constexpr (sizeof(T) > 10) {
#else
#if _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif
        if (sizeof(T) > 10) {
#if _MSC_VER
#pragma warning(pop)
#endif
#endif
          if (isfinite(y)) {
            const T e{(Phi(y) - x)};
            const T u{e * constants<T>::sqrt_2pi * exp(y * y * traits::one_half)};
            y -= u / (1 + y * u * traits::one_half);
          }
        }
        return y;
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_erf(T x) {
        T y{inv_Phi_approx((x + 1) / 2) * constants<T>::one_over_sqrt_2};
        // refinement by Halley rational method
        if (isfinite(y)) {
          const T e{erf(y) - x};
          const T u{e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
          y -= u / (1 + y * u);
        }
        // T is a floating point number type with more than 80 bits, a 2nd iteration is
        // required to reach full machine precision
#if __cplusplus >= 201703L
        if constexpr (sizeof(T) > 10) {
#else
#if _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif
        if (sizeof(T) > 10) {
#if _MSC_VER
#pragma warning(pop)
#endif
#endif
          if (isfinite(y)) {
            const T e{erf(y) - x};
            const T u{e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
            y -= u / (1 + y * u);
          }
        }
        return y;
      }

      template<typename T>
      TRNG_CUDA_ENABLE T inv_erfc(T x) {
        // step size in the Halley step is proportional to erfc, use symmetry to increase
        // numerical accuracy
        const bool flag{x > 1};
        if (flag)
          x = -(x - 1) + 1;
        T y{-inv_Phi_approx(x / 2) * constants<T>::one_over_sqrt_2};
        // refinement by Halley rational method
        if (isfinite(y)) {
          const T e{erfc(y) - x};
          const T u{-e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
          y -= u / (1 + y * u);
        }
        // T is a floating point number type with more than 80 bits, a 2nd iteration is
        // required to reach full machine precision
#if __cplusplus >= 201703L
        if constexpr (sizeof(T) > 10) {
#else
#if _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4127)
#endif
        if (sizeof(T) > 10) {
#if _MSC_VER
#pragma warning(pop)
#endif
#endif
          if (isfinite(y)) {
            const T e{erfc(y) - x};
            const T u{-e * (constants<T>::sqrt_pi_over_2) * exp(y * y)};
            y -= u / (1 + y * u);
          }
        }
        return flag ? -y : y;
      }

    }  // namespace detail

    TRNG_CUDA_ENABLE
    inline float inv_Phi(float x) {
      return detail::inv_Phi<float>(x);
    }

    TRNG_CUDA_ENABLE
    inline double inv_Phi(double x) {
      return detail::inv_Phi<double>(x);
    }

#if !(defined TRNG_CUDA)
    inline long double inv_Phi(long double x) {
      return detail::inv_Phi<long double>(x);
    }
#endif

    // --- inverse of error function  ----------------------------------

    TRNG_CUDA_ENABLE
    inline float inv_erf(float x) {
      return detail::inv_erf(x);
    }

    TRNG_CUDA_ENABLE
    inline double inv_erf(double x) {
      return detail::inv_erf(x);
    }

#if !(defined TRNG_CUDA)
    inline long double inv_erf(long double x) {
      return detail::inv_erf(x);
    }
#endif

    // --- inverse of complementary error function  --------------------

    TRNG_CUDA_ENABLE
    inline float inv_erfc(float x) {
      return detail::inv_erfc(x);
    }

    TRNG_CUDA_ENABLE
    inline double inv_erfc(double x) {
      return detail::inv_erfc(x);
    }

#if !(defined TRNG_CUDA)
    inline long double inv_erfc(long double x) {
      return detail::inv_erfc(x);
    }
#endif

  }  // namespace math

}  // namespace trng

#endif
