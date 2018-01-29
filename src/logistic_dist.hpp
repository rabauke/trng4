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

#if !(defined TRNG_LOGISTIC_DIST_HPP)

#define TRNG_LOGISTIC_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  template<typename float_t=double>
  class logistic_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type theta_, eta_;
    public:
      TRNG_CUDA_ENABLE
      result_type theta() const { return theta_; }
      TRNG_CUDA_ENABLE
      void theta(result_type theta_new) { theta_=theta_new; }
      TRNG_CUDA_ENABLE
      result_type eta() const { return eta_; }
      TRNG_CUDA_ENABLE
      void eta(result_type eta_new) { eta_=eta_new; }
      TRNG_CUDA_ENABLE
      param_type() : theta_(1), eta_(0) {
      }
      TRNG_CUDA_ENABLE
      param_type(result_type theta, result_type eta) : theta_(theta), eta_(eta) {
      }

      friend class logistic_dist;

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
            << p.theta() << ' ' << p.eta() 
            << ')';
        out.flags(flags);
        return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
                 param_type &p) {
        float_t theta, eta;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed |
                 std::ios_base::left);
        in >> utility::delim('(')
           >> theta >> utility::delim(' ')
           >> eta >> utility::delim(')');
        if (in)
          p=param_type(theta, eta);
        in.flags(flags);
        return in;
      }
      
    };
    
  private:
    param_type p;
    
    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      return p.eta()-math::ln((1-x)/x)*p.theta();
    }
    
  public:
    // constructor
    TRNG_CUDA_ENABLE
    logistic_dist(result_type theta, result_type eta) : p(theta, eta) {
    }
    TRNG_CUDA_ENABLE
    explicit logistic_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() { }
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r) {
      return icdf_(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE
    result_type operator()(R &r, const param_type &p) {
      logistic_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return -math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p=p_new; }
    TRNG_CUDA_ENABLE
    result_type theta() const { return p.theta(); }
    TRNG_CUDA_ENABLE
    void theta(result_type theta_new) { p.theta(theta_new); }
    TRNG_CUDA_ENABLE
    result_type eta() const { return p.eta(); }
    TRNG_CUDA_ENABLE
    void eta(result_type eta_new) { p.eta(eta_new); }
    // probability density function  
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      result_type t1(math::exp(-math::abs((x-p.eta())/p.theta())));
      result_type t2(1+t1);
      return t1/(p.theta()*t2*t2);
    }
    // cumulative density function 
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      result_type t1(math::exp(-math::abs((x-p.eta())/p.theta())));
      return x>=p.eta() ? (1/(1+t1)) : (t1/(1+t1));
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
      if (x==0)
        return -math::numeric_limits<result_type>::infinity();
      if (x==1)
        return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const typename logistic_dist<float_t>::param_type &p1, 
			 const typename logistic_dist<float_t>::param_type &p2) {
    return p1.theta()==p2.theta() and p1.eta()==p2.eta();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const typename logistic_dist<float_t>::param_type &p1, 
			 const typename logistic_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }  
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator==(const logistic_dist<float_t> &g1, 
			 const logistic_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE
  inline bool operator!=(const logistic_dist<float_t> &g1, 
			 const logistic_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const logistic_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[logistic " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     logistic_dist<float_t> &g) {
    typename logistic_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[logistic ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
