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

#if !(defined TRNG_GAMMA_DIST_HPP)

#define TRNG_GAMMA_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  template<typename float_t=double>
  class gamma_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type kappa_, theta_;
    public:
      TRNG_CUDA_ENABLE
      result_type kappa() const { return kappa_; }
      TRNG_CUDA_ENABLE
      void kappa(result_type kappa_new) { kappa_=kappa_new; }
      TRNG_CUDA_ENABLE
      result_type theta() const { return theta_; }
      TRNG_CUDA_ENABLE
      void theta(result_type theta_new) { theta_=theta_new; }
      TRNG_CUDA_ENABLE
      param_type() : kappa_(1), theta_(1) {
      }
      TRNG_CUDA_ENABLE
      explicit param_type(result_type kappa, result_type theta) : kappa_(kappa), theta_(theta) {
      }

      friend class gamma_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
		 const param_type &p) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed |
		  std::ios_base::left);
	out << '('
	    << std::setprecision(math::numeric_limits<result_type>::digits10+1)	  
	    << p.kappa() << ' ' << p.theta() 
	    << ')';
	out.flags(flags);
	return out;
      }
  
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t kappa, theta;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> kappa >> utility::delim(' ')
	   >> theta >> utility::delim(')');
	if (in)
	  p=param_type(kappa, theta);
	in.flags(flags);
	return in;
      }
    };
    
  private:
    param_type p;

    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      if (x<=math::numeric_limits<result_type>::epsilon())
	return 0;
      if (p.kappa()==1)  // special case of exponential distribution
	return -math::ln(1-x)*p.theta();
      result_type ln_Gamma_kappa=math::ln_Gamma(p.kappa());
      result_type y=p.kappa(), y_old;
      int num_iterations=0;
      do {
	++num_iterations;
	y_old=y;
	result_type f0=math::GammaP(p.kappa(), y)-x;
	result_type f1=math::exp((p.kappa()-1)*math::ln(y)-y-ln_Gamma_kappa);
	result_type f2=f1*(p.kappa()-1-y)/y;
	y-=f0/f1*(1+f0*f2/(2*f1*f1));
      } while (num_iterations<16 &&
	       math::abs((y-y_old)/y)>16*math::numeric_limits<result_type>::epsilon());
      return y*p.theta();
    }
    
  public:
    // constructor
    TRNG_CUDA_ENABLE
    gamma_dist(result_type kappa, result_type theta) : p(kappa, theta) {
    }
    TRNG_CUDA_ENABLE
    explicit gamma_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() { }
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r) {
      return icdf_(utility::uniformco<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r, const param_type &p) {
      gamma_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return 0; }
    TRNG_CUDA_ENABLE
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p=p_new; }
    TRNG_CUDA_ENABLE
    result_type kappa() const { return p.kappa(); }
    TRNG_CUDA_ENABLE
    void kappa(result_type kappa_new) { p.kappa(kappa_new); }
    TRNG_CUDA_ENABLE
    result_type theta() const { return p.theta(); }
    TRNG_CUDA_ENABLE
    void theta(result_type theta_new) { p.theta(theta_new); }
    // probability density function  
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      if (x<0)
	return 0;
      else {
	x/=p.theta();
	return math::exp((p.kappa()-1)*math::ln(x)-x-math::ln_Gamma(p.kappa())) /(p.theta());
      }
    }
    // cumulative density function 
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      if (x<=0)
	return 0;
      else
	return math::GammaP(p.kappa(), x/p.theta());
    }
    // inverse cumulative density function 
    TRNG_CUDA_ENABLE
    result_type icdf(result_type x) const {
      if (x<=0 or x>=1) {
#if !(defined __CUDA_ARCH__)
	errno=EDOM;
#endif
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      if (x==0)
	return 0;
      if (x==1)
	return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const typename gamma_dist<float_t>::param_type &p1, 
			 const typename gamma_dist<float_t>::param_type &p2) {
    return p1.kappa()==p2.kappa() and p1.theta()==p2.theta();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const typename gamma_dist<float_t>::param_type &p1, 
			 const typename gamma_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const gamma_dist<float_t> &g1, 
			 const gamma_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const gamma_dist<float_t> &g1, 
			 const gamma_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const gamma_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[gamma " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     gamma_dist<float_t> &g) {
    typename gamma_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[gamma ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
