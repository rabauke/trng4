// Copyright (C) 2007 Heiko Bauke <heiko.bauke@physik.uni-magdeburg.de>
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

#if !(defined TRNG_RAYLEIGH_DIST_HPP)

#define TRNG_RAYLEIGH_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class rayleigh_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double nu_;
    public:
      double nu() const { return nu_; }
      void nu(double nu_new) { nu_=nu_new; }
      explicit param_type(double nu) : nu_(nu) {
      }
      friend class rayleigh_dist;
    };
    
  private:
    param_type p;
   
  public:
    // constructor
    explicit rayleigh_dist(double nu) : p(nu) {
    }
    explicit rayleigh_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return icdf(utility::uniformoo(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      rayleigh_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0.0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double nu() const { return p.nu(); }
    void nu(double nu_new) { p.nu(nu_new); }
    // probability density function  
    double pdf(double x) const {
      if (x<=0.0)
	return 0.0;
      double t=x/(p.nu()*p.nu());
      return t*math::exp(-0.5*t*x);
    }
    // cumulative density function 
    double cdf(double x) const {
      if (x<=0.0)
	return 0.0;
      return 1.0-math::exp(-0.5*x*x/(p.nu()*p.nu()));
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      return p.nu()*math::sqrt(-2.0*math::ln(1.0-x));
    }
  };
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const rayleigh_dist::param_type &p1, 
			 const rayleigh_dist::param_type &p2) {
    return p1.nu()==p2.nu();
  }
  inline bool operator!=(const rayleigh_dist::param_type &p1, 
			 const rayleigh_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const rayleigh_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.nu()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     rayleigh_dist::param_type &p) {
    double nu;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> nu >> utility::delim(')');
    if (in)
      p=rayleigh_dist::param_type(nu);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const rayleigh_dist &g1, 
			 const rayleigh_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const rayleigh_dist &g1, 
			 const rayleigh_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const rayleigh_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[rayleigh " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     rayleigh_dist &g) {
    rayleigh_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[rayleigh ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
