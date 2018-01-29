// Copyright (c) 2000-2018, Heiko Bauke
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

#if !(defined TRNG_TENT_DIST_HPP)

#define TRNG_TENT_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // non-uniform random number generator class
  template<typename float_t=double>
  class tent_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type m_, d_;
    public:
      TRNG_CUDA_ENABLE
      result_type m() const { return m_; }
      TRNG_CUDA_ENABLE
      void m(result_type m_new) { m_=m_new; }
      TRNG_CUDA_ENABLE
      result_type d() const { return d_; }
      TRNG_CUDA_ENABLE
      void d(result_type d_new) { d_=d_new; }
      TRNG_CUDA_ENABLE
      param_type() : m_(0), d_(1) {
      }
      TRNG_CUDA_ENABLE
      param_type(result_type m, result_type d) : m_(m), d_(d) {
      }

      friend class tent_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
		 const param_type &p) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed |
		  std::ios_base::left);
	out << '('
	    << std::setprecision(math::numeric_limits<float_t>::digits10+1) 
	    << p.m() << ' ' << p.d()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t m, d;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> m >> utility::delim(' ')
	   >> d >> utility::delim(')');
	if (in)
	  p=param_type(m, d);
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      return x<result_type(1)/result_type(2) ?
        (-1+math::sqrt(    2*x))*p.d()+p.m() :
        ( 1-math::sqrt(2-2*x))*p.d()+p.m();
    }

  public:
    // constructor
    TRNG_CUDA_ENABLE
    tent_dist(result_type m, result_type d) : p(m, d) {
    }
    TRNG_CUDA_ENABLE
    explicit tent_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() { }
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r) {
      return icdf_(utility::uniformcc<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r, const param_type &p) {
      tent_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return p.m()-p.d(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return p.m()+p.d(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p=p_new; }
    TRNG_CUDA_ENABLE
    result_type m() const { return p.m(); }
    TRNG_CUDA_ENABLE
    void m(result_type m_new) { p.m(m_new); }
    TRNG_CUDA_ENABLE
    result_type d() const { return p.d(); }
    TRNG_CUDA_ENABLE
    void d(result_type d_new) { p.d(d_new); }
    // probability density function  
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      x-=p.m();
      if (x<=-p.d() or x>=p.d())
        return 0;
      return x<0 ? (1+x/p.d())/p.d() : (1-x/p.d())/p.d();
    }
    // cumulative density function 
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      x-=p.m();
      if (x<=-p.d())
        return 0;
      if (x<=0)
        return (x+p.d())*(x+p.d())/(result_type(2)*p.d()*p.d());
      if (x<p.d())
        return 1-(x-p.d())*(x-p.d())/(result_type(2)*p.d()*p.d());
      return 1;
    }
    // inverse cumulative density function 
    TRNG_CUDA_ENABLE
    result_type icdf(result_type x) const {
      if (x<0 or x>1) {
#if !(defined __CUDA_ARCH__)
	errno=EDOM;
#endif
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      return icdf_(x);
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const typename tent_dist<float_t>::param_type &p1, 
			 const typename tent_dist<float_t>::param_type &p2) {
    return p1.m()==p2.m() and p1.d()==p2.d();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const typename tent_dist<float_t>::param_type &p1, 
			 const typename tent_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const tent_dist<float_t> &g1, 
			 const tent_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const tent_dist<float_t> &g1, 
			 const tent_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const tent_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[tent " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     tent_dist<float_t> &g) {
    typename tent_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[tent ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
