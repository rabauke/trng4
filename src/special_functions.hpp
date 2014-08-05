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

#if !(defined TRNG_SPECIAL_FUNCTIONS_HPP)

#define TRNG_SPECIAL_FUNCTIONS_HPP

#include <trng/config.hpp>
#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/math.hpp>
#include <trng/constants.hpp>
#include <trng/utility.hpp>
#include <cerrno>
#include <algorithm>

namespace trng {

  namespace math {

    // --- log-Gamma function ------------------------------------------

    namespace detail {
    
      template<typename T>
      struct ln_Gamma_traits;
    
      template<>
      struct ln_Gamma_traits<float> {
	static float one_half() throw() { 
	  return 0.5f;
	}
	static float ln_sqrt_2pi() throw() {
	  return 0.91893853320467274177f;
	}
	static float b(int i) throw() {
	  // B_{2k}/(2k(2k-1))   k=1, 2, ...
	  static const float b_[]={1.0f/12.0f,
				   -1.0f/360.0f,
				   1.0f/1260.0f};
	  return b_[i];
	}
	static const int size_b=3;
      };
    
      template<>
      struct ln_Gamma_traits<double> {
	static double one_half() throw() { 
	  return 0.5;
	}
	static double ln_sqrt_2pi() throw() {
	  return 0.91893853320467274177;
	}
	static double b(int i) throw() {
	  // B_{2k}/(2k(2k-1))   k=1, 2, ...
	  static const double b_[]={1.0/12.0,
				    -1.0/360.0,
				    1.0/1260.0,
				    -1.0/1680.0,
				    1.0/1188.0,
				    -691.0/360360.0,
				    1.0/156.0};
	  return b_[i];
	}
	static const int size_b=7;
      };

      template<>
      struct ln_Gamma_traits<long double> {
	static long double one_half() throw() { 
	  return 0.5l;
	}
	static long double ln_sqrt_2pi() throw() {
	  return 0.91893853320467274177l;
	}
	static long double b(int i) throw() {
	  // B_{2k}/(2k(2k-1))   k=1, 2, ...
	  static const long double b_[]={1.0l/12.0l,
					 -1.0l/360.0l,
					 1.0l/1260.0l,
					 -1.0l/1680.0l,
					 1.0l/1188.0l,
					 -691.0l/360360.0l,
					 1.0l/156.0l,
					 -3617.0l/122400.0l,
					 43867.0l/125400.0l};
	  return b_[i];
	}
	static const int size_b=9;
      };

      // calculate ln(Gamma(x+1)) for x near zero using Legendre's series
      template<typename T>
      T ln_Gamma_1(T x) {
	static const T eps(T(4)*numeric_limits<T>::epsilon());
	if (abs(x)<eps)
	  return T(0);
	T t(x*constants<T>::pi());
	T sum( ln_Gamma_traits<T>::one_half()*
	       ln( t*(T(1)-x)/(sin(t)*(T(1)+x)) )  );
	sum+=ln_Gamma_traits<T>::one_min_gamma()*x;  // (1-gamma)*x
	T x2(x*x);
	int i(0);
	while (i<ln_Gamma_traits<T>::size_a) {
	  x*=x2;
	  t=ln_Gamma_traits<T>::a(i++)*x;
	  if (abs(t)<eps*sum)
	    break;
	  sum+=t;
	}
	return sum;
      }

      // calculate ln(Gamma(x)) by Lanczos' approximation
      template<typename T>
      T ln_Gamma_Lanczos(T x);

      template<>
      inline float ln_Gamma_Lanczos<float>(float x) {
	return ln(+0.299890266072888E-2f/(x+4.0f)
		  -0.308748865044984E1f/(x+3.0f)
		  +0.6019440944395479E2f/(x+2.0f)
		  -0.2168366808191931E3f/(x+1.0f)
		  +0.190955171863804E3f/x
		  +0.250662827022856E1f)
	  -(x+4.5f)+(x+0.5f)*ln(x+4.5f);
      }
    
      template<>
      inline double ln_Gamma_Lanczos<double>(double x) {
	return ln(-0.1710538478644311E-5/(x+7.0)
		  +0.8683645856906762E-1/(x+6.0)
		  -0.1567563009175129E2/(x+5.0)
		  +0.3873696975776843E3/(x+4.0)
		  -0.2900131673187631E4/(x+3.0)
		  +0.8679533396416264E4/(x+2.0)
		  -0.1102260304013762E5/(x+1.0)
		  +0.4951528076618453E4/x
		  +0.2506628274630859E1)
	  -(x+7.5)+(x-0.5)*ln(x+7.5);
      }
    
      template<>
      inline long double ln_Gamma_Lanczos<long double>(long double x) {
	return ln(-0.1008550193581785E-7l/(x+9.0l)
		  +0.1349602861619936E-6l/(x+8.0l)
		  -0.1860783595745135E-1l/(x+7.0l)
		  +0.6531512536073623E1l/(x+6.0l)
		  -0.2711587882903386E3l/(x+5.0l)
		  +0.3262648132327785E4l/(x+4.0l)
		  -0.1591247789342777E5l/(x+3.0l)
		  +0.3582345988044706E5l/(x+2.0l)
		  -0.3713646057463839E5l/(x+1.0l)
		  +0.1432889034103444E5l/x
		  +0.2506628274631001E1l)
	  -(x+8.5l)+(x-0.5l)*ln(x+8.5l);
      }

      // calculate ln(Gamma(x)) for large x using an asympotic series
      template<typename T>
      T ln_Gamma_infty(T x) {
	static const T eps(T(4)*numeric_limits<T>::epsilon());
	T sum((x-ln_Gamma_traits<T>::one_half())*ln(x)-x+
	      ln_Gamma_traits<T>::ln_sqrt_2pi());
	x=T(1)/x;
	T x2(x*x), t;
	int i(0);
	while (i<ln_Gamma_traits<T>::size_b) {
	  t=ln_Gamma_traits<T>::b(i++)*x;
	  if (abs(t)<eps*sum)
	    break;
	  x*=x2;
	  sum+=t;
	}
	return sum;
      }
    
      // calculate ln(Gamma(x)) for positive x
      template<typename T>
      inline T ln_Gamma(T x) {
	if (x<T(20)) 
	  return ln_Gamma_Lanczos(x);
	return ln_Gamma_infty(x);
      }
    
    }

    TRNG_CUDA_ENABLE
    inline float ln_Gamma(float x) {
#if defined(TRNG_HAVE_LGAMMAF) or defined(TRNG_CUDA)
      return ::lgammaf(x);
#elif __cplusplus >= 201103L
      return std::lgamma(x);
#else
      if (x>0.0f)
	return detail::ln_Gamma(x);
      float t(abs(sin(constants<float>::pi()*x)));
      if (t<4.0f*numeric_limits<float>::epsilon())
	return numeric_limits<float>::infinity();
      return ln(constants<float>::pi()/(-x*t))-detail::ln_Gamma(-x);
#endif
    }
  
    TRNG_CUDA_ENABLE
    inline double ln_Gamma(double x) {
#if defined(TRNG_HAVE_LGAMMA) or defined(TRNG_CUDA)
      return ::lgamma(x);
#elif __cplusplus >= 201103L
      return std::lgamma(x);
#else
      if (x>0.0)
	return detail::ln_Gamma(x);
      double t(abs(sin(constants<double>::pi()*x)));
      if (t<4.0*numeric_limits<double>::epsilon())
	return numeric_limits<double>::infinity();
      return ln(constants<double>::pi()/(-x*t))-detail::ln_Gamma(-x);
#endif
    }

#if !(defined __CUDA_ARCH__)
    inline long double ln_Gamma(long double x) {
#if defined(TRNG_HAVE_LGAMMAL)
      return ::lgammal(x);
#elif __cplusplus >= 201103L
      return std::lgamma(x);
#else
      if (x>0.0l)
	return detail::ln_Gamma(x);
      long double t(abs(sin(constants<long double>::pi()*x)));
      if (t<4.0l*numeric_limits<long double>::epsilon())
	return numeric_limits<long double>::infinity();
      return ln(constants<long double>::pi()/(-x*t))-detail::ln_Gamma(-x);
#endif
    }
#endif

    // --- Gamma function ----------------------------------------------

    namespace detail {

      template<typename T>
      struct Gamma_traits;
    
      template<>
      struct Gamma_traits<float> {
	static float one_half() throw() { 
	  return 0.5f;
	}
	static float lim() throw() { 
	  return 20.0f;
	}
	static float sqrt_2pi() throw() {
	  return 2.506628274631000502416f;
	}
	static float a(int i) throw() {
	  static const float a_[]={1.0f,
				   1.0f/12.0f,
				   1.0f/288.0f};
	  return a_[i];
	}
	static const int size_a=3;
      };
    
      template<>
      struct Gamma_traits<double> {
	static double one_half() throw() { 
	  return 0.5;
	}
	static double lim() throw() { 
	  return 20.0;
	}
	static double sqrt_2pi() throw() {
	  return 2.506628274631000502416;
	}
	static double a(int i) throw() {
	  static const double a_[]={1.0,
				    1.0/12.0,
				    1.0/288.0,
				    -139.0/51840.0,
				    -571.0/2488320.0,
				    163879.0/209018880.0,
				    5246819.0/75246796800.0,
				    -534703531.0/902961561600.0,
				    -4483131259.0/86684309913600.0};
	  return a_[i];
	}
	static const int size_a=9;
      };
    
      template<>
      struct Gamma_traits<long double> {
	static long double one_half() throw() { 
	  return 0.5l;
	}
	static long double lim() throw() { 
	  return 20.0l;
	}
	static long double sqrt_2pi() throw() {
	  return 2.506628274631000502416l;
	}
	static long double a(int i) throw() {
	  static const long double a_[]={1.0l,
					 1.0l/12.0l,
					 1.0l/288.0l,
					 -139.0l/51840.0l,
					 -571.0l/2488320.0l,
					 163879.0l/209018880.0l,
					 5246819.0l/75246796800.0l,
					 -534703531.0l/902961561600.0l,
					 -4483131259.0l/86684309913600.0l,
					 432261921612371.0l/514904800886784000.0l,
					 6232523202521089.0l/86504006548979712000.0l};
	  return a_[i];
	}
	static const int size_a=11;
      };

      // calculate Gamma(x) by Lanczos' approximation
      template<typename T>
      T Gamma_Lanczos(T x);

      template<>
      inline float Gamma_Lanczos<float>(float x) {
	return (+0.299890266072888E-2f/(x+4.0f)
		-0.308748865044984E1f/(x+3.0f)
		+0.6019440944395479E2f/(x+2.0f)
		-0.2168366808191931E3f/(x+1.0f)
		+0.190955171863804E3f/x
		+0.250662827022856E1f)*
	  exp(-(x+4.5f)+(x+0.5f)*ln(x+4.5f));
      }

      template<>
      inline double Gamma_Lanczos<double>(double x) {
	return (-0.1710538478644311E-5/(x+7.0)
		+0.8683645856906762E-1/(x+6.0)
		-0.1567563009175129E2/(x+5.0)
		+0.3873696975776843E3/(x+4.0)
		-0.2900131673187631E4/(x+3.0)
		+0.8679533396416264E4/(x+2.0)
		-0.1102260304013762E5/(x+1.0)
		+0.4951528076618453E4/x
		+0.2506628274630859E1)*
	  exp(-(x+7.5)+(x-0.5)*ln(x+7.5));
      }
    
      template<>
      inline long double Gamma_Lanczos<long double>(long double x) {
	return (-0.1008550193581785E-7l/(x+9.0l)
		+0.1349602861619936E-6l/(x+8.0l)
		-0.1860783595745135E-1l/(x+7.0l)
		+0.6531512536073623E1l/(x+6.0l)
		-0.2711587882903386E3l/(x+5.0l)
		+0.3262648132327785E4l/(x+4.0l)
		-0.1591247789342777E5l/(x+3.0l)
		+0.3582345988044706E5l/(x+2.0l)
		-0.3713646057463839E5l/(x+1.0l)
		+0.1432889034103444E5l/x
		+0.2506628274631001E1l)*
	  exp(-(x+8.5l)+(x-0.5l)*ln(x+8.5l));
      }

      // calculate Gamma(x) for positive x
      template<typename T>
      T Gamma(T x) {
	static const T eps(T(4)*numeric_limits<T>::epsilon());
	if (x<Gamma_traits<T>::lim())
	  return Gamma_Lanczos(x);
	// use Stirling series
	T x1(T(1)/x), x2(T(1)), sum(0), t;
	int i(0);
	while (i<Gamma_traits<T>::size_a) {
	  t=Gamma_traits<T>::a(i++)*x2;
	  if (abs(t)<eps*sum)
	    break;
	  x2*=x1;
	  sum+=t;
	}
	return sum*Gamma_traits<T>::sqrt_2pi()*
	  pow(x1, Gamma_traits<T>::one_half()-x)*exp(-x);
      }

    }

    TRNG_CUDA_ENABLE
    inline float Gamma(float x) {
#if defined(TRNG_HAVE_TGAMMAF) or defined(TRNG_CUDA)
      return ::tgammaf(x);
#elif __cplusplus >= 201103L
      return std::tgamma(x);
#else
      if (x>0.0f)
	return detail::Gamma(x);
      float t(sin(x*constants<float>::pi()));
      if (abs(t)<4.0f*numeric_limits<float>::epsilon())
	return (t>=0.0f ? 1.0f : -1.0f)*numeric_limits<float>::infinity();
      return constants<float>::pi()/(-x*detail::Gamma(-x)*t);
#endif
    }

    TRNG_CUDA_ENABLE
    inline double Gamma(double x) {
#if defined(TRNG_HAVE_TGAMMA) or defined(TRNG_CUDA)
      return ::tgamma(x);
#elif __cplusplus >= 201103L
      return std::tgamma(x);
#else
      if (x>0.0)
	return detail::Gamma(x);
      double t(sin(x*constants<double>::pi()));
      if (abs(t)<4.0*numeric_limits<double>::epsilon())
	return (t>=0.0 ? 1.0 : -1.0)*numeric_limits<double>::infinity();
      return constants<double>::pi()/(-x*detail::Gamma(-x)*t);
#endif
    }

#if !(defined __CUDA_ARCH__)
    inline long double Gamma(long double x) {
#if defined(TRNG_HAVE_TGAMMAL)
      return ::tgammal(x);
#elif __cplusplus >= 201103L
      return std::tgamma(x);
#else
      if (x>0.0l)
	return detail::Gamma(x);
      long double t(sin(x*constants<long double>::pi()));
      if (abs(t)<4.0l*numeric_limits<long double>::epsilon())
	return (t>=0.0l ? 1.0l : -1.0l)*numeric_limits<long double>::infinity();
      return constants<long double>::pi()/(-x*detail::Gamma(-x)*t);
#endif
    }
#endif

    // --- Beta function -----------------------------------------------

    TRNG_CUDA_ENABLE
    inline float Beta(float x, float y) {
      return exp(ln_Gamma(x)+ln_Gamma(y)-ln_Gamma(x+y));
    }

    TRNG_CUDA_ENABLE
    inline double Beta(double x, double y) {
      return exp(ln_Gamma(x)+ln_Gamma(y)-ln_Gamma(x+y));
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double Beta(long double x, long double y) {
      return exp(ln_Gamma(x)+ln_Gamma(y)-ln_Gamma(x+y));
    }
#endif

    // --- ln of binomial coefficient ----------------------------------

    TRNG_CUDA_ENABLE
    inline float ln_binomial(float n, float m) {
      return ln_Gamma(n+1.f)-ln_Gamma(m+1.f)-ln_Gamma(n-m+1.f);
    }

    TRNG_CUDA_ENABLE
    inline double ln_binomial(double n, double m) {
      return ln_Gamma(n+1.)-ln_Gamma(m+1.)-ln_Gamma(n-m+1.);
    }

#if !(defined __CUDA_ARCH__)
    inline long double ln_binomial(long double n, long double m) {
      return ln_Gamma(n+1.l)-ln_Gamma(m+1.l)-ln_Gamma(n-m+1.l);
    }
#endif

    // --- Pochhammer function -----------------------------------------

    inline float Pochhammer(float x, float a) {
      return exp(ln_Gamma(x+a)-ln_Gamma(x));
    }

    inline double Pochhammer(double x, double a) {
      return exp(ln_Gamma(x+a)-ln_Gamma(x));
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double Pochhammer(long double x, long double a) {
      return exp(ln_Gamma(x+a)-ln_Gamma(x));
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
      template<typename T>
      TRNG_CUDA_ENABLE
      T GammaP_ser(T a, T x, bool by_Gamma_a) {
	const int itmax=32;
	 const T eps=T(4)*numeric_limits<T>::epsilon();
	if (x<eps)
	  return T(0);
	T xx(T(1)/a), n(a), sum(xx);
	int i(0);
	do {
	  ++n;
	  ++i;
	  xx*=x/n;
	  sum+=xx;
	} while (abs(xx)>eps*abs(sum) and i<itmax);
	if (by_Gamma_a)
	  return exp(-x+a*ln(x)-math::ln_Gamma(a))*sum;
	return exp(-x+a*ln(x))*sum;
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
      template<typename T>
      TRNG_CUDA_ENABLE
      T GammaQ_cf(T a, T x, bool by_Gamma_a) {
	const T itmax=T(32);
	const T eps=T(4)*numeric_limits<T>::epsilon();
	const T min=T(4)*numeric_limits<T>::min();
	// set up for evaluating continued fraction by modied Lentz's method
	T del, ai, bi(x+T(1)-a), ci(T(1)/min), di(T(1)/bi), h(di), i(T(0));
	do { // iterate 
	  ++i;
	  ai=-i*(i-a);
	  bi+=T(2);
	  di=ai*di+bi;
	  if (abs(di)<min) 
	    di=min; 
	  ci=bi+ai/ci; 
	  if (abs(ci)<min)
	    ci=min;
	  di=T(1)/di;
	  del=di*ci;
	  h*=del;
	} while ((abs(del-T(1))>eps) and i<itmax);
	if (by_Gamma_a)
	  return exp(-x+a*ln(x)-math::ln_Gamma(a))*h;
	return exp(-x+a*ln(x))*h;
      }
    
      // P(a, x) and gamma(a, x)
      template<typename T>
      TRNG_CUDA_ENABLE
      T GammaP(T a, T x, bool by_Gamma_a) {
	if (x<T(0) or a<=T(0))
	  return numeric_limits<T>::signaling_NaN();
	if (by_Gamma_a) {
	  if (x<a+T(1))
	    return GammaP_ser(a, x, true);
	  return T(1)-GammaQ_cf(a, x, true);
	}
	if (x<a+T(1))
	  return GammaP_ser(a, x, false);
	return math::Gamma(a)-GammaQ_cf(a, x, false);
      }

      // Q(a, x) and Gamma(a, x)
      template<typename T>
      TRNG_CUDA_ENABLE
      T GammaQ(T a, T x, bool by_Gamma_a) {
	if (x<T(0) or a<=T(0))
	  return numeric_limits<T>::signaling_NaN();
	if (by_Gamma_a) {
	  if (x<a+T(1))
	    return T(1)-GammaP_ser(a, x, true);
	  return GammaQ_cf(a, x, true);
	}
	if (x<a+T(1))
	  return math::Gamma(a)-GammaP_ser(a, x, false);
	return GammaQ_cf(a, x, false);
      }
    
    }
  
    // P(x, a)
    TRNG_CUDA_ENABLE
    inline float GammaP(float a, float x) {
      return detail::GammaP(a, x, true);
    }

    TRNG_CUDA_ENABLE
    inline double GammaP(double a, double x) {
      return detail::GammaP(a, x, true);
    }

#if !(defined __CUDA_ARCH__)
    inline long double GammaP(long double a, long double x) {
      return detail::GammaP(a, x, true);
    }
#endif

    // Q(x, a)
    TRNG_CUDA_ENABLE
    inline float GammaQ(float a, float x) {
      return detail::GammaQ(a, x, true);
    }

    TRNG_CUDA_ENABLE
    inline double GammaQ(double a, double x) {
      return detail::GammaQ(a, x, true);
    }

#if !(defined __CUDA_ARCH__)
    inline long double GammaQ(long double a, long double x) {
      return detail::GammaQ(a, x, true);
    }
#endif

    // gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float inc_gamma(float a, float x) {
      return detail::GammaP(a, x, false);
    }

    TRNG_CUDA_ENABLE
    inline double inc_gamma(double a, double x) {
      return detail::GammaP(a, x, false);
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double inc_gamma(long double a, long double x) {
      return detail::GammaP(a, x, false);
    }
#endif

    // Gamma(x, a)
    TRNG_CUDA_ENABLE
    inline float cinc_gamma(float a, float x) {
      return detail::GammaQ(a, x, false);
    }
  
    TRNG_CUDA_ENABLE
    inline double cinc_gamma(double a, double x) {
      return detail::GammaQ(a, x, false);
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double cinc_gamma(long double a, long double x) {
      return detail::GammaQ(a, x, false);
    }
#endif

    namespace detail {
      
      template<typename T>
      TRNG_CUDA_ENABLE
      T inv_GammaP(T a, T p) {
        const T eps=sqrt(numeric_limits<T>::epsilon()),
          a1=a-T(1), glna=math::ln_Gamma(a), 
          lna1=ln(a1), afac=exp(a1*(lna1-T(1))-glna);
        T x, t;
        // initial guess
        if (a>T(1)) {
          const T pp=p<T(1)/T(2) ? p : T(1)-p;
          t=sqrt(-T(2)*ln(pp));
          x=static_cast<T>((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481))-t);
          x=p<T(1)/T(2) ? -x : x;
          x=utility::max(T(1)/T(1000),
		     a*pow(T(1)-T(1)/(T(9)*a)-x/(T(3)*sqrt(a)), T(3)));
        } else {
          t=static_cast<T>(1.0-a*(0.253+a*0.12));
          x=p<t ? (pow(p/t, T(1)/a)) : (T(1)-ln(T(1)-(p-t)/(T(1)-t)));
        }
        // refinement by Halley's method
        for (int i=0; i<16; ++i) {
          if (x<T(0)) {
            x=T(0);
            break;
          }
          const T err=GammaP(a, x, true)-p;
          if (a>T(1)) 
            t=afac*exp(-(x-a1)+a1*(ln(x)-lna1));
          else
            t=exp(-x+a1*ln(x)-glna);
          const T u=err/t;
          t=u/(T(1)-utility::min(T(1), u*((a-T(1))/x-T(1)))/T(2));
          x-=t;
          x=x<=T(0) ? (x+t)/T(2) : x;
          if (abs(t)<eps*x)
            break;
        }
        return x;
      }
      
    }

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
      TRNG_CUDA_ENABLE
      T Beta_I(T x, T p, T q, T norm) {
        if (p<=0 or q<=0 or x<0 or x>1) {
#if !(defined __CUDA_ARCH__)
          errno=EDOM;
#endif
          return numeric_limits<T>::quiet_NaN();
        }
        const T eps=4*numeric_limits<T>::epsilon();
        T psq=p+q, cx=1-x;
        bool flag=(p<psq*x);
        if (flag) {
          // use  I(x, p, q) = 1-I(1-x, q, p)
	  utility::swap(x, cx);
	  utility::swap(p, q);
        }
        T term=1, i=1, y=1, rx=x/cx, temp=q-i;
        int s=static_cast<int>(q+cx*psq);
        if (s==0)
          rx=x;
        while (true) {
          term*=temp*rx/(p+i);
          y+=term;
          temp=abs(term);
          if (temp<=eps and temp<=eps*y)
            break;
          i++;
          s--;
          if (s>=0) {
            temp=q-i;
            if (s==0)
              rx=x;
          } else {
            temp=psq;
            psq++;
          }
        }
        y*=exp(p*ln(x)+(q-1)*ln(cx))/p/norm;
        if (flag)
          y=1-y;
        return y;
      }
      
    }
    
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

//     // see Applied Statistics (1973), vol.22, no.3, pp.411--414
//     // algorithm AS 64
//     namespace detail {

//       TRNG_CUDA_ENABLE
//       inline float inv_Beta_I_approx(float x, float p, float q, float norm) {
//         float r=sqrt(-2.0f*ln(x));
//         float y=r-(2.30753f+0.27061f*r)/(1.0f+(0.99229f+0.04481f*r)*r);
//         r=q+q;
//         float t=1.0f/(9.0f*q);
//         t=r*(1.0f-t+y*sqrt(t)); 
//         t=t*t*t;
//         if (t>0.0f) {
//           t=(4.0f*p+r-2.0f)/t;
//           if (t>1.0f)
//             y=1.0f-2.0f/(t+1.0f);
//           else
//             y=exp(ln(x*p*norm)/p);
//         } else
//           y=1.0f-exp(ln((1.0f-x)*q*norm)/q);
//         return y;
//       }

//       TRNG_CUDA_ENABLE
//       inline double inv_Beta_I_approx(double x, double p, double q, double norm) {
//         double r=sqrt(-2.0*ln(x));
//         double y=r-(2.30753+0.27061*r)/(1.0+(0.99229+0.04481*r)*r);
//         r=q+q;
//         double t=1.0/(9.0*q);
//         t=r*(1.0-t+y*sqrt(t)); 
//         t=t*t*t;
//         if (t>0.0) {
//           t=(4.0*p+r-2.0)/t;
//           if (t>1.0)
//             y=1.0-2.0/(t+1.0);
//           else
//             y=exp(ln(x*p*norm)/p);
//         } else
//           y=1.0-exp(ln((1.0-x)*q*norm)/q);
//         return y;
//       }

//       inline long double inv_Beta_I_approx(long double x, long double p, long double q, long double norm) {
//         long double r=sqrt(-2.0l*ln(x));
//         long double y=r-(2.30753l+0.27061l*r)/(1.0l+(0.99229l+0.04481l*r)*r);
//         r=q+q;
//         long double t=1.0l/(9.0l*q);
//         t=r*(1.0l-t+y*sqrt(t)); 
//         t=t*t*t;
//         if (t>0.0l) {
//           t=(4.0l*p+r-2.0l)/t;
//           if (t>1.0l)
//             y=1.0l-2.0l/(t+1.0l);
//           else
//             y=exp(ln(x*p*norm)/p);
//         } else
//           y=1.0l-exp(ln((1.0l-x)*q*norm)/q);
//         return y;
//       }

//       template<typename T>
//       TRNG_CUDA_ENABLE
//       T inv_Beta_I(T x, T p, T q, T norm) {
//         if (p<=0 or q<=0 or x<0 or x>1) {
// #if !(defined __CUDA_ARCH__)
//           errno=EDOM;
// #endif
//           return numeric_limits<T>::quiet_NaN();
//         }
//         const T eps=4*numeric_limits<T>::epsilon();
//         bool flag;
//         if (flag=(2*x<1)) {
//           // use  I(x, p, q) = 1-I(1-x, q, p)
//           x=1-x;
// 	  utility::swap(p, q);
//         }
//         // initial approximation
//         T y=inv_Beta_I_approx(x, p, q, norm);     
//         // Newton-Rapson iteration
//         T r=1-p, t=1-q, yy;
// 	int iters=0;
// 	T err=T(2), err_old;
//         do {
//           yy=Beta_I(y, p, q, norm);
//           yy=(yy-x)*norm*exp(r*ln(y)+t*ln(1-y));
//           y-=yy;
// 	  err_old=err;
// 	  err=abs(yy);
// 	  ++iters;
//         } while (err>eps and err_old>err and iters<32);
//         if (flag)
//           y=1-y;
//         return y;
//       }
//     }

    namespace detail {

      template<typename T>
      TRNG_CUDA_ENABLE
      inline T inv_Beta_I(T x, T p, T q, T norm) {
	// solve via Householder method
	T y(p/(q+p-1));
	for (int i=0; i<math::numeric_limits<T>::digits; ++i) {
	  T f(math::Beta_I(y, p, q, norm)-x);
	  T df(math::pow(1-y, q-1)*math::pow(y, p-1));
	  df/=norm;
	  T ddf(-math::pow(1-y, q-1)*(q-1)/(1-y)*math::pow(y, p)+
		math::pow(1-y, q-1)*math::pow(y, p)*p/y);
	  ddf/=norm;
	  T dy(f/df*(1+f*ddf/(2*df*df)));
	  y-=dy;
	  if (math::abs(dy)<math::numeric_limits<T>::epsilon())
	    break;
	}
	return y;
      }
      
    }

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
    inline long double inv_Beta_I(long double x, long double p, long double q, long double norm) {
      return detail::inv_Beta_I(x, p, q, norm);
    }

    inline long double inv_Beta_I(long double x, long double p, long double q) {
      return detail::inv_Beta_I(x, p, q, Beta(p, q));
    }
#endif

    // --- error function ----------------------------------------------

    TRNG_CUDA_ENABLE
    inline float erf(float x) {
#if defined(TRNG_HAVE_ERFF) or defined(TRNG_CUDA)
      return ::erff(x);
#elif __cplusplus >= 201103L
      return std::erf(x);
#else
      return x<0.0f ? -GammaP(0.5f, x*x) : GammaP(0.5f, x*x);
#endif
    }

    TRNG_CUDA_ENABLE
    inline double erf(double x) {
#if defined(TRNG_HAVE_ERF) or defined(TRNG_CUDA)
      return ::erf(x);
#elif __cplusplus >= 201103L
      return std::erf(x);
#else
      return x<0.0 ? -GammaP(0.5, x*x) : GammaP(0.5, x*x);
#endif
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double erf(long double x) {
#if defined(TRNG_HAVE_ERFL)
      return ::erfl(x);
#elif __cplusplus >= 201103L
      return std::erf(x);
#else
      return x<0.0l ? -GammaP(0.5l, x*x) : GammaP(0.5l, x*x);
#endif
    }
#endif
    
    // --- complementary error function --------------------------------

    TRNG_CUDA_ENABLE
    inline float erfc(float x) {
#if defined(TRNG_HAVE_ERFCF) or defined(TRNG_CUDA)
      return ::erfcf(x);
#elif __cplusplus >= 201103L
      return std::erfc(x);
#else
      return x<0.0f ? 1.0f+GammaP(0.5f, x*x) : GammaQ(0.5f, x*x);
#endif
    }

    TRNG_CUDA_ENABLE
    inline double erfc(double x) {
#if defined(TRNG_HAVE_ERFC) or defined(TRNG_CUDA)
      return ::erfc(x);
#elif __cplusplus >= 201103L
      return std::erfc(x);
#else
      return x<0.0 ? 1.0+GammaP(0.5, x*x) : GammaQ(0.5, x*x);
#endif
    }
  
#if !(defined __CUDA_ARCH__)
    inline long double erfc(long double x) {
#if defined(TRNG_HAVE_ERFCL)
      return ::erfcl(x);
#elif __cplusplus >= 201103L
      return std::erfc(x);
#else
      return x<0.0l ? 1.0l+GammaP(0.5l, x*x) : GammaQ(0.5l, x*x);
#endif
    }
#endif

    // --- normal distribution function  -------------------------------

    TRNG_CUDA_ENABLE
    inline float Phi(float x) {
      return 0.5f+0.5f*erf(constants<float>::one_over_sqrt_2()*x);
    }

    TRNG_CUDA_ENABLE
    inline double Phi(double x) {
      return 0.5+0.5*erf(constants<double>::one_over_sqrt_2()*x);
    }

    inline long double Phi(long double x) {
      return 0.5l+0.5l*erf(constants<long double>::one_over_sqrt_2()*x);
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
	static float a(int i) throw() {
	  const float a_[]={
	    -3.969683028665376e+01f,   2.209460984245205e+02f,
	    -2.759285104469687e+02f,   1.383577518672690e+02f,
	    -3.066479806614716e+01f,   2.506628277459239e+00f};
	  return a_[i];
	}
	TRNG_CUDA_ENABLE
	static float b(int i) throw() {
	  const float b_[]={
	    -5.447609879822406e+01f,   1.615858368580409e+02f,
	    -1.556989798598866e+02f,   6.680131188771972e+01f,
	    -1.328068155288572e+01f};
	  return b_[i];
	}
	TRNG_CUDA_ENABLE
	static float c(int i) throw() {
	  const float c_[]={
	    -7.784894002430293e-03f,  -3.223964580411365e-01f,
	    -2.400758277161838e+00f,  -2.549732539343734e+00f,
	    4.374664141464968e+00f,    2.938163982698783e+00f};
	  return c_[i];
	}
	TRNG_CUDA_ENABLE
	static float d(int i) throw() {
	  const float d_[]={
	    7.784695709041462e-03f,    3.224671290700398e-01f,
	    2.445134137142996e+00f,    3.754408661907416e+00f};
	  return d_[i];
	}
	TRNG_CUDA_ENABLE
	static float x_low() throw() {
	  return 0.02425f;
	}
	TRNG_CUDA_ENABLE
	static float x_high() throw() {
	  return 1.0f-0.02425f;
	}
	TRNG_CUDA_ENABLE
	static float zero() throw() {
	  return 0.0f;
	}
	TRNG_CUDA_ENABLE
	static float one() throw() {
	  return 1.0f;
	}
	TRNG_CUDA_ENABLE
	static float one_half() throw() {
	  return 0.5f;
	}
	TRNG_CUDA_ENABLE
	static float minus_two() throw() {
	  return -2.0f;
	}
      };
      
      template<>
      struct inv_Phi_traits<double> {
	TRNG_CUDA_ENABLE
	static double a(int i) throw() {
	  const double a_[]={
	    -3.969683028665376e+01,   2.209460984245205e+02,
	    -2.759285104469687e+02,   1.383577518672690e+02,
	    -3.066479806614716e+01,   2.506628277459239e+00};
	  return a_[i];
	}
	TRNG_CUDA_ENABLE
	static double b(int i) throw() {
	  const double b_[]={
	    -5.447609879822406e+01,   1.615858368580409e+02,
	    -1.556989798598866e+02,   6.680131188771972e+01,
	    -1.328068155288572e+01};
	  return b_[i];
	}
	TRNG_CUDA_ENABLE
	static double c(int i) throw() {
	  const double c_[]={
	    -7.784894002430293e-03,  -3.223964580411365e-01,
	    -2.400758277161838e+00,  -2.549732539343734e+00,
	    4.374664141464968e+00,    2.938163982698783e+00};
	  return c_[i];
	}
	TRNG_CUDA_ENABLE
	static double d(int i) throw() {
	  const double d_[]={
	    7.784695709041462e-03,    3.224671290700398e-01,
	    2.445134137142996e+00,    3.754408661907416e+00};
	  return d_[i];
	}
	TRNG_CUDA_ENABLE
	static double x_low() throw() {
	  return 0.02425;
	}
	TRNG_CUDA_ENABLE
	static double x_high() throw() {
	  return 1.0-0.02425;
	}
	TRNG_CUDA_ENABLE
	static double zero() throw() {
	  return 0.0;
	}
	TRNG_CUDA_ENABLE
	static double one() throw() {
	  return 1.0;
	}
	TRNG_CUDA_ENABLE
	static double one_half() throw() {
	  return 0.5;
	}
	TRNG_CUDA_ENABLE
	static double minus_two() throw() {
	  return -2.0;
	}
      };

      template<>
      struct inv_Phi_traits<long double> {
	TRNG_CUDA_ENABLE
	static long double a(int i) throw() {
	  const long double a_[]={
	    -3.969683028665376e+01l,   2.209460984245205e+02l,
	    -2.759285104469687e+02l,   1.383577518672690e+02l,
	    -3.066479806614716e+01l,   2.506628277459239e+00l};
	  return a_[i];
	}
	TRNG_CUDA_ENABLE
	static long double b(int i) throw() {
	  const long double b_[]={
	    -5.447609879822406e+01l,   1.615858368580409e+02l,
	    -1.556989798598866e+02l,   6.680131188771972e+01l,
	    -1.328068155288572e+01l};
	  return b_[i];
	}
	TRNG_CUDA_ENABLE
	static long double c(int i) throw() {
	  const long double c_[]={
	    -7.784894002430293e-03l,  -3.223964580411365e-01l,
	    -2.400758277161838e+00l,  -2.549732539343734e+00l,
	    4.374664141464968e+00l,    2.938163982698783e+00l};
	  return c_[i];
	}
	TRNG_CUDA_ENABLE
	static long double d(int i) throw() {
	  const long double d_[]={
	    7.784695709041462e-03l,    3.224671290700398e-01l,
	    2.445134137142996e+00l,    3.754408661907416e+00l};
	  return d_[i];
	}
	TRNG_CUDA_ENABLE
	static long double x_low() throw() {
	  return 0.02425l;
	}
	TRNG_CUDA_ENABLE
	static long double x_high() throw() {
	  return 1.0l-0.02425l;
	}
	TRNG_CUDA_ENABLE
	static long double zero() throw() {
	  return 0.0l;
	}
	TRNG_CUDA_ENABLE
	static long double one() throw() {
	  return 1.0l;
	}
	TRNG_CUDA_ENABLE
	static long double one_half() throw() {
	  return 0.5l;
	}
	TRNG_CUDA_ENABLE
	static long double minus_two() throw() {
	  return -2.0l;
	}
      };

      template<typename T>
      TRNG_CUDA_ENABLE
      T inv_Phi(T x) {
	if (x<inv_Phi_traits<T>::zero() or x>inv_Phi_traits<T>::one()) {
#if !(defined __CUDA_ARCH__)
	  errno=EDOM;
#endif
	  return numeric_limits<T>::quiet_NaN();
	} 
	if (x==inv_Phi_traits<T>::zero())
	  return -numeric_limits<T>::infinity();
	if (x==inv_Phi_traits<T>::one())
	  return numeric_limits<T>::infinity();
	T t, q;
	if (x<inv_Phi_traits<T>::x_low()) {
	  // Rational approximation for lower region
	  q=sqrt(inv_Phi_traits<T>::minus_two()*ln(x));
	  t=(((((inv_Phi_traits<T>::c(0)*q + inv_Phi_traits<T>::c(1))*q +
		inv_Phi_traits<T>::c(2))*q + inv_Phi_traits<T>::c(3))*q +
	      inv_Phi_traits<T>::c(4))*q + inv_Phi_traits<T>::c(5)) /
	    ((((inv_Phi_traits<T>::d(0)*q + inv_Phi_traits<T>::d(1))*q +
	       inv_Phi_traits<T>::d(2))*q + inv_Phi_traits<T>::d(3))*q +
	     inv_Phi_traits<T>::one());
	} else if (x<inv_Phi_traits<T>::x_high()) {
	  // Rational approximation for central region
	  q=x-inv_Phi_traits<T>::one_half();
	  T r=q*q;
	  t=(((((inv_Phi_traits<T>::a(0)*r + inv_Phi_traits<T>::a(1))*r + 
		inv_Phi_traits<T>::a(2))*r + inv_Phi_traits<T>::a(3))*r + 
	      inv_Phi_traits<T>::a(4))*r + inv_Phi_traits<T>::a(5))*q /
	    (((((inv_Phi_traits<T>::b(0)*r + inv_Phi_traits<T>::b(1))*r +
		inv_Phi_traits<T>::b(2))*r + inv_Phi_traits<T>::b(3))*r +
	      inv_Phi_traits<T>::b(4))*r + inv_Phi_traits<T>::one());
	} else {
	  // Rational approximation for upper region
	  q=sqrt(inv_Phi_traits<T>::minus_two()*ln(1.0-x));
	  t=-(((((inv_Phi_traits<T>::c(0)*q + inv_Phi_traits<T>::c(1))*q +
		 inv_Phi_traits<T>::c(2))*q + inv_Phi_traits<T>::c(3))*q +
	       inv_Phi_traits<T>::c(4))*q + inv_Phi_traits<T>::c(5)) /
	    ((((inv_Phi_traits<T>::d(0)*q + inv_Phi_traits<T>::d(1))*q +
	       inv_Phi_traits<T>::d(2))*q + inv_Phi_traits<T>::d(3))*q +
	     inv_Phi_traits<T>::one());	
	}
	// refinement by Halley rational method
	if (numeric_limits<T>::epsilon()<1e-9) {
	  T e(Phi(t)-x);
	  T u(e*constants<T>::sqrt_2pi()*exp(t*t*inv_Phi_traits<T>::one_half()));
	  t-=u/(inv_Phi_traits<T>::one()+t*u*inv_Phi_traits<T>::one_half());
	}
	return t;
      }
      
    }

    TRNG_CUDA_ENABLE
    inline float inv_Phi(float x) {
      return detail::inv_Phi<float>(x);
    }
    
    TRNG_CUDA_ENABLE
    inline double inv_Phi(double x) {
      return detail::inv_Phi<double>(x);
    }

#if !(defined __CUDA_ARCH__)
    inline long double inv_Phi(long double x) {
      return detail::inv_Phi<long double>(x);
    }
#endif

    // --- inverse of error function  ----------------------------------
    
    // see http://mathworld.wolfram.com/InverseErf.html
    // see The On-Line Encyclopedia of Integer Sequences! 
    // http://www.research.att.com/~njas/sequences/A007019
    // http://www.research.att.com/~njas/sequences/A092676

    TRNG_CUDA_ENABLE
    inline float inv_erf(float x) {
      if (abs(x)<1.0f/8.0f) {
	x*=0.886226925452758013649085f;  // sqrt(pi)/2
	float x2=x*x, x3=x2*x, x4=x2*x2;
	return x + (1.0f/3.0f + 7.0f/30.0f*x2 + 127.0f/630.0f*x4)*x3;
      }
      return inv_Phi(0.5f*(x+1.0f))*constants<float>::one_over_sqrt_2();
    }
    
    TRNG_CUDA_ENABLE
    inline double inv_erf(double x) {
      if (abs(x)<1.0/20.0) {
	x*=0.886226925452758013649085;  // sqrt(pi)/2
	double x2=x*x, x3=x2*x, x4=x2*x2, x5=x3*x2;
	return x + (1.0/3.0 + 127.0/630.0*x4)*x3 +
	  (7.0/30.0 + 4369.0/22680.0*x4)*x5;
      }
      return inv_Phi(0.5*(x+1.0))*constants<double>::one_over_sqrt_2();
    }
    
#if !(defined __CUDA_ARCH__)
    inline long double inv_erf(long double x) {
      if (abs(x)<1.0l/24.0l) {
	x*=0.886226925452758013649085l;  // sqrt(pi)/2
	long double x2=x*x, x3=x2*x, x4=x2*x2, x5=x3*x2, x7=x3*x4;
	return x + 1.0l/3.0l*x3 + 7.0l/30.0l*x5 + 
	  127.0l/630.0l*x7 + 4369.0l/22680.0l*x4*x5 + 
	  34807.0l/178200.0l*x4*x7;
      }
      return inv_Phi(0.5l*(x+1.0l))*constants<long double>::one_over_sqrt_2();
    }
#endif

    // --- inverse of complementary error function  --------------------
    
    TRNG_CUDA_ENABLE
    inline float inv_erfc(float x) {
      return -inv_Phi(0.5f*x)*constants<float>::one_over_sqrt_2();
    }
    
    TRNG_CUDA_ENABLE
    inline double inv_erfc(double x) {
      return -inv_Phi(0.5*x)*constants<double>::one_over_sqrt_2();
    }
    
#if !(defined __CUDA_ARCH__)
    inline long double inv_erfc(long double x) {
      return -inv_Phi(0.5l*x)*constants<long double>::one_over_sqrt_2();
    }
#endif
    
  }

}

#endif
