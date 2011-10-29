// Copyright (C) 2000-2007 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#if !(defined TRNG_EXTREME_VALUE_DIST_HPP)

#define TRNG_EXTREME_VALUE_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class extreme_value_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double theta_, eta_;
    public:
      double theta() const { return theta_; }
      void theta(double theta_new) { theta_=theta_new; }
      double eta() const { return eta_; }
      void eta(double eta_new) { eta_=eta_new; }
      explicit param_type(double theta, double eta) : theta_(theta), eta_(eta) {
      }
      friend class extreme_value_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit extreme_value_dist(double theta, double eta) : p(theta, eta) {
    }
    explicit extreme_value_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return p.eta()+p.theta()*math::ln(-math::ln(utility::uniformoo(r)));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      extreme_value_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return -math::numeric_limits<double>::infinity(); }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double theta() const { return p.theta(); }
    void theta(double theta_new) { p.theta(theta_new); }
    double eta() const { return p.eta(); }
    void eta(double eta_new) { p.eta(eta_new); }
    // probability density function  
    double pdf(double x) const {
      x-=p.eta();
      x/=p.theta();
      return math::exp(x-math::exp(x))/p.theta();
    }
    // cumulative density function 
    double cdf(double x) const {
      x-=p.eta();
      x/=p.theta();
      return 1.0-math::exp(-math::exp(x));
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<=0.0 || x>=1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      if (x==0.0)
	return -math::numeric_limits<double>::infinity();
      if (x==1.0)
	return math::numeric_limits<double>::infinity();
      return p.eta()+p.theta()*math::ln(-math::ln(1.0-x));
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const extreme_value_dist::param_type &p1, 
			 const extreme_value_dist::param_type &p2) {
    return p1.theta()==p2.theta() && p1.eta()==p2.eta();
  }
  inline bool operator!=(const extreme_value_dist::param_type &p1, 
			 const extreme_value_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const extreme_value_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.theta() << ' ' << p.eta() 
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     extreme_value_dist::param_type &p) {
    double theta, eta;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> theta >> utility::delim(' ')
       >> eta >> utility::delim(')');
    if (in)
      p=extreme_value_dist::param_type(theta, eta);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const extreme_value_dist &g1, 
			 const extreme_value_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const extreme_value_dist &g1, 
			 const extreme_value_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const extreme_value_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[extreme_value " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     extreme_value_dist &g) {
    extreme_value_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[extreme_value ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
