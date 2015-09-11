// Copyright (c) 2000-2015, Heiko Bauke
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

#if !(defined TRNG_NEGATIVE_BINOMIAL_DIST_HPP)

#define TRNG_NEGATIVE_BINOMIAL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>

namespace trng {

  // non-uniform random number generator class
  class negative_binomial_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      double p_;
      int r_;
      std::vector<double> P_;
      
      void calc_probabilities() {
      	P_=std::vector<double>();
       	int x=0;
       	double p=0.0;
       	while (p<1.0-1.0/4096.0) {
       	  p+=math::exp(math::ln_Gamma(static_cast<double>(r_+x))-
		       math::ln_Gamma(static_cast<double>(x+1))-
		       math::ln_Gamma(static_cast<double>(r_)))*
	    math::pow(p_, r_)*math::pow(1-p_, x);
	  P_.push_back(p);
	  std::cerr << x << '\t' << p << '\n';
	  ++x;
      	}
	P_.push_back(1);
      }
      
    public:
      double p() const { return p_; }
      void p(double p_new) { p_=p_new;  calc_probabilities(); }
      int r() const { return r_; }
      void r(int r_new) { r_=r_new;  calc_probabilities(); }
      param_type() :
        p_(0), r_(0) {
      }
      param_type(double p, int r) :
        p_(p), r_(r) {
        calc_probabilities();
      }
      friend class negative_binomial_dist;
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    explicit negative_binomial_dist(double p, int r) : P(p, r) {
    }
    explicit negative_binomial_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      double p(utility::uniformco<double>(r));
      int x(utility::discrete(p, P.P_.begin(), P.P_.end()));
      if (x+1==P.P_.size()) {
	p-=cdf(x);
        while (p>0) {
          ++x;
	  p-=pdf(x);
	}
      }
      return x;
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      negative_binomial_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return P.r(); }
    int max() const { return math::numeric_limits<int>::max(); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    double p() const { return P.p(); }
    void p(double p_new) { P.p(p_new); }
    double r() const { return P.r(); }
    void r(double r_new) { P.r(r_new); }
    // probability density function  
    double pdf(int x) const {
      return x<0 ? 0.0 : math::exp(math::ln_Gamma(static_cast<double>(P.r()+x))-
				   math::ln_Gamma(static_cast<double>(x+1))-
				   math::ln_Gamma(static_cast<double>(P.r())))*
	       math::pow(P.p(), P.r())*math::pow(1-P.p(), x);
    }
    // cumulative density function 
    double cdf(int x) const {
      double res=0;
      while (x>=0) {
	res+=pdf(x);
	--x;
      }
      return res;
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const negative_binomial_dist::param_type &p1, 
			 const negative_binomial_dist::param_type &p2) {
    return p1.p()==p2.p() and p1.r()==p2.r();
  }
  inline bool operator!=(const negative_binomial_dist::param_type &p1, 
			 const negative_binomial_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const negative_binomial_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << P.p() << ' ' << P.r()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     negative_binomial_dist::param_type &P) {
    double p;
    int r;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> p >> utility::delim(' ')
       >> r >> utility::delim(')');
    if (in)
      P=negative_binomial_dist::param_type(p, r);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const negative_binomial_dist &g1, 
			 const negative_binomial_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const negative_binomial_dist &g1, 
			 const negative_binomial_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const negative_binomial_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[negative_binomial " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     negative_binomial_dist &g) {
    negative_binomial_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[negative_binomial ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif

