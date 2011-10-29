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

#if !(defined TRNG_STUDENT_T_DIST_HPP)

#define TRNG_STUDENT_T_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class student_t_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      int nu_;
    public:
      int nu() const { return nu_; }
      void nu(int nu_new) { nu_=nu_new; }
      explicit param_type(int nu) : nu_(nu) {
      }
      friend class student_t_dist;
    };
    
  private:
    param_type p;

    // inverse cumulative density function
    double icdf_(double x) const {
      double t=math::inv_Beta_I(x, 0.5*p.nu(), 0.5*p.nu());
      return math::sqrt(p.nu()/(t*(1.0-t)))*(t-0.5);
    }

  public:
    // constructor
    explicit student_t_dist(int nu) : p(nu) {
    }
    explicit student_t_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return icdf_(utility::uniformoo(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      student_t_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return -math::numeric_limits<double>::infinity(); }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    int nu() const { return p.nu(); }
    void nu(int nu_new) { p.nu(nu_new); }
    // probability density function  
    double pdf(double x) const {
      double norm=
	math::exp(math::ln_Gamma(0.5*(p.nu()+1))-
		  math::ln_Gamma(0.5*p.nu()))/
	math::sqrt(math::constants<double>::pi()*p.nu());
      return norm*math::pow(1.0+x*x/p.nu(), -0.5*(p.nu()+1));
    }
    // cumulative density function 
    double cdf(double x) const {
      double t1=+math::sqrt(x*x+p.nu());
      double t2=0.5*(x+t1)/t1;
      return math::Beta_I(t2, 0.5*p.nu(), 0.5*p.nu());
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
      return icdf_(x);
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const student_t_dist::param_type &p1, 
			 const student_t_dist::param_type &p2) {
    return p1.nu()==p2.nu();
  }
  inline bool operator!=(const student_t_dist::param_type &p1, 
			 const student_t_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const student_t_dist::param_type &p) {
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
	     student_t_dist::param_type &p) {
    int nu;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> nu >> utility::delim(')');
    if (in)
      p=student_t_dist::param_type(nu);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const student_t_dist &g1, 
			 const student_t_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const student_t_dist &g1, 
			 const student_t_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const student_t_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[student_t " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     student_t_dist &g) {
    student_t_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[student_t ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
