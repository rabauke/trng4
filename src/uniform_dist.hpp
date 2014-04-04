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

#if !(defined TRNG_UNIFORM_DIST_HPP)

#define TRNG_UNIFORM_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  template<typename float_t=double>
  class uniform_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type a_, b_, d_;
      TRNG_CUDA_ENABLE
      result_type d() const { return d_; }
    public:
      TRNG_CUDA_ENABLE
      result_type a() const { return a_; }
      TRNG_CUDA_ENABLE
      void a(result_type a_new) { a_=a_new; d_=b_-a_; }
      TRNG_CUDA_ENABLE
      result_type b() const { return b_; }
      TRNG_CUDA_ENABLE
      void b(result_type b_new) { b_=b_new; d_=b_-a_; }
      TRNG_CUDA_ENABLE
      param_type() :
	a_(0), b_(1), d_(1) {
      }
      TRNG_CUDA_ENABLE
      param_type(result_type a, result_type b) :
	a_(a), b_(b), d_(b-a) {
      }

      friend class uniform_dist<float_t>;

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
	    << p.a() << ' ' << p.b()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t a, b;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> a >> utility::delim(' ')
	   >> b >> utility::delim(')');
	if (in)
	  p=param_type(a, b);
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    TRNG_CUDA_ENABLE
    uniform_dist(result_type a, result_type b) : p(a, b) {
    }
    TRNG_CUDA_ENABLE
    explicit uniform_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() { }
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r) {
      return p.d()*utility::uniformco<result_type>(r)+p.a();
    }
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r, const param_type &p) {
      uniform_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return p.a(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return p.b(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p=p_new; }
    TRNG_CUDA_ENABLE
    result_type a() const { return p.a(); }
    TRNG_CUDA_ENABLE
    void a(result_type a_new) { p.a(a_new); }
    TRNG_CUDA_ENABLE
    result_type b() const { return p.b(); }
    TRNG_CUDA_ENABLE
    void b(result_type b_new) { p.b(b_new); }
    // probability density function  
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      if (x<p.a() or x>=p.b())
	return 0.0;
      return 1.0/p.d();
    }
    // cumulative density function 
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      if (x<p.a())
	return 0;
      if (x>=p.b())
	return 1;
      return (x-p.a())/p.d();
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
      return x*p.d()+p.a(); 
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const typename uniform_dist<float_t>::param_type &p1, 
			 const typename uniform_dist<float_t>::param_type &p2) {
    return p1.a()==p2.a() and p1.b()==p2.b();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const typename uniform_dist<float_t>::param_type &p1, 
			 const typename uniform_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const uniform_dist<float_t> &g1, 
			 const uniform_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const uniform_dist<float_t> &g1, 
			 const uniform_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // -------------------------------------------------------------------

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[uniform " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform_dist<float_t> &g) {
    typename uniform_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[uniform ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
