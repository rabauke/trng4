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

#if !(defined TRNG_GEOMETRIC_DIST_HPP)

#define TRNG_GEOMETRIC_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <climits>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>

namespace trng {

  // non-uniform random number generator class
  class geometric_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      double p_, q_, one_over_ln_q_;

      TRNG_CUDA_ENABLE
      double q() const {  return q_;  }
      TRNG_CUDA_ENABLE
      double one_over_ln_q() const {  return one_over_ln_q_;  }

    public:
      TRNG_CUDA_ENABLE
      double p() const { return p_; }
      TRNG_CUDA_ENABLE
      void p(double p_new) { 
	p_=p_new;  q_=1.0-p_;  one_over_ln_q_=1.0/math::ln(q_);
      }
      TRNG_CUDA_ENABLE
      explicit param_type(double p=0.5) :
	p_(p), q_(1.0-p_), one_over_ln_q_(1.0/math::ln(q_)) {
      }
      friend class geometric_dist;
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    TRNG_CUDA_ENABLE
    explicit geometric_dist(double p) : P(p) {
    }
    TRNG_CUDA_ENABLE
    explicit geometric_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() { }
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE
    int operator()(R &r) {
      return static_cast<int>(math::ln(utility::uniformoo<double>(r))*
			      P.one_over_ln_q());
    }
    template<typename R>
    TRNG_CUDA_ENABLE
    int operator()(R &r, const param_type &p) {
      geometric_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    int min() const { return 0; }
    TRNG_CUDA_ENABLE
    int max() const { return INT_MAX; }
    TRNG_CUDA_ENABLE
    param_type param() const { return P; }
    TRNG_CUDA_ENABLE
    void param(const param_type &P_new) { P=P_new; }
    TRNG_CUDA_ENABLE
    double p() const { return P.p(); }
    TRNG_CUDA_ENABLE
    void p(double p_new) { P.p(p_new); }
    // probability density function  
    TRNG_CUDA_ENABLE
    double pdf(int x) const {
      return x<0 ? 0.0 : P.p()*math::pow(P.q(), static_cast<double>(x));
    }
    // cumulative density function 
    TRNG_CUDA_ENABLE
    double cdf(int x) const {
      return x<0 ? 0.0 : 1.0-math::pow(P.q(), static_cast<double>(x+1));
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  TRNG_CUDA_ENABLE
  inline bool operator==(const geometric_dist::param_type &p1, 
			 const geometric_dist::param_type &p2) {
    return p1.p()==p2.p();
  }

  TRNG_CUDA_ENABLE
  inline bool operator!=(const geometric_dist::param_type &p1, 
			 const geometric_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const geometric_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << P.p()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     geometric_dist::param_type &P) {
    double p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> p >> utility::delim(')');
    if (in)
      P=geometric_dist::param_type(p);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  TRNG_CUDA_ENABLE
  inline bool operator==(const geometric_dist &g1, 
			 const geometric_dist &g2) {
    return g1.param()==g2.param();
  }

  TRNG_CUDA_ENABLE
  inline bool operator!=(const geometric_dist &g1, 
			 const geometric_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const geometric_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[geometric " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     geometric_dist &g) {
    geometric_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[geometric ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif

