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

#if !(defined TRNG_CAUCHY_DIST_HPP)

#define TRNG_CAUCHY_DIST_HPP

#include <trng/constants.hpp>
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
  class cauchy_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type theta_, eta_;
    public:
      result_type theta() const { return theta_; }
      void theta(result_type theta_new) { theta_=theta_new; }
      result_type eta() const { return eta_; }
      void eta(result_type eta_new) { eta_=eta_new; }
      param_type() : theta_(1), eta_(0) {
      }
      param_type(result_type theta, result_type eta) : theta_(theta), eta_(eta) {
      }

      friend class cauchy_dist;

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
    result_type icdf_(result_type x) const {
      result_type t=x*math::constants<result_type>::pi();
      return p.eta()-p.theta()*math::cos(t)/math::sin(t);
    }
    
  public:
    // constructor
    cauchy_dist(result_type theta, result_type eta) : p(theta, eta) {
    }
    explicit cauchy_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return icdf_(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      cauchy_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return -math::numeric_limits<result_type>::infinity(); }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type theta() const { return p.theta(); }
    void theta(result_type theta_new) { p.theta(theta_new); }
    result_type eta() const { return p.eta(); }
    void eta(result_type eta_new) { p.eta(eta_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      x-=p.eta();
      x/=p.theta();
      return math::constants<result_type>::one_over_pi()/(1.0+x*x)/p.theta();
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      x-=p.eta();
      x/=p.theta();
      return math::constants<result_type>::one_over_pi()*math::atan(x)+
	result_type(1)/result_type(2);
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<=0 or x>=1) {
	errno=EDOM;
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
  inline bool operator==(const typename cauchy_dist<float_t>::param_type &p1, 
			 const typename cauchy_dist<float_t>::param_type &p2) {
    return p1.theta()==p2.theta() and p1.eta()==p2.eta();
  }

  template<typename float_t>
  inline bool operator!=(const typename cauchy_dist<float_t>::param_type &p1, 
			 const typename cauchy_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const cauchy_dist<float_t> &g1, 
			 const cauchy_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }
  
  template<typename float_t>
  inline bool operator!=(const cauchy_dist<float_t> &g1, 
			 const cauchy_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const cauchy_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[cauchy " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     cauchy_dist<float_t> &g) {
    typename cauchy_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[cauchy ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
