// Copyright (C) 2000-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
//  
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License in
// version 2 as published by the Free Software Foundation.
//  
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
//  

#if !(defined TRNG_WEIBULL_DIST_HPP)

#define TRNG_WEIBULL_DIST_HPP

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
  class weibull_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type theta_, beta_;
    public:
      result_type theta() const { return theta_; }
      void theta(result_type theta_new) { theta_=theta_new; }
      result_type beta() const { return beta_; }
      void beta(result_type beta_new) { beta_=beta_new; }
      param_type() : theta_(1), beta_(1) {
      }
      param_type(result_type theta, result_type beta) : theta_(theta), beta_(beta) {
      }

      friend class weibull_dist;

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
	    << p.theta() << ' ' << p.beta() 
	    << ')';
	out.flags(flags);
	return out;
      }
  
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t theta, beta;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> theta >> utility::delim(' ')
	   >> beta >> utility::delim(')');
	if (in)
	  p=param_type(theta, beta);
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    weibull_dist(result_type theta, result_type beta) : p(theta, beta) {
    }
    explicit weibull_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return p.theta()*
	math::pow(-math::ln(utility::uniformoc<result_type>(r)), 1/p.beta());
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      weibull_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return 0; }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type theta() const { return p.theta(); }
    void theta(result_type theta_new) { p.theta(theta_new); }
    result_type beta() const { return p.beta(); }
    void beta(result_type beta_new) { p.beta(beta_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      if (x<0)
        return 0;
      x/=p.theta();
      if (x>0) {
        result_type t(math::pow(x, p.beta()));
        return p.beta()*t/x*math::exp(-t);
      }
      // x==0
      if (p.beta()==1)
        return 1/p.theta();
      if (p.beta()>1)
        return 0;
      return math::numeric_limits<result_type>::quiet_NaN();
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      x/=p.theta();
      if (x<=0)
        return 0;
      return 1-math::exp(-math::pow(x, p.beta()));
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<0 or x>=1) {
        errno=EDOM;
        return math::numeric_limits<result_type>::quiet_NaN();
      }
      return p.theta()*math::pow(-math::ln(1-x), 1/p.beta());
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename weibull_dist<float_t>::param_type &p1, 
			 const typename weibull_dist<float_t>::param_type &p2) {
    return p1.theta()==p2.theta() and p1.beta()==p2.beta();
  }

  template<typename float_t>
  inline bool operator!=(const typename weibull_dist<float_t>::param_type &p1, 
			 const typename weibull_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const weibull_dist<float_t> &g1, 
			 const weibull_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const weibull_dist<float_t> &g1, 
			 const weibull_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const weibull_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[weibull " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     weibull_dist<float_t> &g) {
    typename weibull_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[weibull ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
