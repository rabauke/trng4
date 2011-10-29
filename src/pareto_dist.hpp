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

#if !(defined TRNG_PARETO_DIST_HPP)

#define TRNG_PARETO_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class pareto_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double gamma_, theta_;
    public:
      double gamma() const { return gamma_; }
      void gamma(double gamma_new) { gamma_=gamma_new; }
      double theta() const { return theta_; }
      void theta(double theta_new) { theta_=theta_new; }
      explicit param_type(double gamma, double theta) : gamma_(gamma), theta_(theta) {
      }
      friend class pareto_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit pareto_dist(double gamma, double theta) : p(gamma, theta) {
    }
    explicit pareto_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return p.theta()*math::pow(utility::uniformoc(r), -1.0/p.gamma())-
	p.theta();
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      pareto_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double gamma() const { return p.gamma(); }
    void gamma(double gamma_new) { p.gamma(gamma_new); }
    double theta() const { return p.theta(); }
    void theta(double theta_new) { p.theta(theta_new); }
    // probability density function  
    double pdf(double x) const {
      if (x<0.0)
	return 0.0;
      else
	return p.gamma()/p.theta()*
	  math::pow(1+x/p.theta(), -p.gamma()-1.0);
    }
    // cumulative density function 
    double cdf(double x) const {
      if (x<=0.0)
	return 0.0;
      else
	return 1.0-math::pow(1+x/p.theta(), -p.gamma());
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<=0.0 || x>=1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      if (x==0.0)
	return 0.0;
      if (x==1.0)
	return math::numeric_limits<double>::infinity();
      return p.theta()*math::pow(1.0-x, -p.gamma()-1.0)+p.theta();
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const pareto_dist::param_type &p1, 
			 const pareto_dist::param_type &p2) {
    return p1.gamma()==p2.gamma() && p1.theta()==p2.theta();
  }
  inline bool operator!=(const pareto_dist::param_type &p1, 
			 const pareto_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const pareto_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.gamma() << ' ' << p.theta() 
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     pareto_dist::param_type &p) {
    double gamma, theta;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> gamma >> utility::delim(' ')
       >> theta >> utility::delim(')');
    if (in)
      p=pareto_dist::param_type(gamma, theta);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const pareto_dist &g1, 
			 const pareto_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const pareto_dist &g1, 
			 const pareto_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const pareto_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[pareto " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     pareto_dist &g) {
    pareto_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[pareto ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
