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

#if !(defined TRNG_TENT_DIST_HPP)

#define TRNG_TENT_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // non-uniform random number generator class
  class tent_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double m_, d_;
    public:
      double m() const { return m_; }
      void m(double m_new) { m_=m_new; }
      double d() const { return d_; }
      void d(double d_new) { d_=d_new; }
      explicit param_type(double m, double d) :	m_(m), d_(d) {
      }
      friend class tent_dist;
    };
    
  private:
    param_type p;

    double icdf_(double x) const {
      return x<0.5 ?
        (-1.0+math::sqrt(    2.0*x))*p.d()+p.m() :
        ( 1.0-math::sqrt(2.0-2.0*x))*p.d()+p.m();
    }

  public:
    // constructor
    explicit tent_dist(double m, double d) : p(m, d) {
    }
    explicit tent_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return icdf_(utility::uniformcc(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      tent_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return p.m()-p.d(); }
    double max() const { return p.m()+p.d(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double m() const { return p.m(); }
    void m(double m_new) { p.m(m_new); }
    double d() const { return p.d(); }
    void d(double d_new) { p.d(d_new); }
    // probability density function  
    double pdf(double x) const {
      x-=p.m();
      if (x<=-p.d() || x>=p.d())
        return 0.0;
      return x<0.0 ? (1.0+x/p.d())/p.d() : (1.0-x/p.d())/p.d();
    }
    // cumulative density function 
    double cdf(double x) const {
      x-=p.m();
      if (x<=p.d())
        return 0;
      if (x<=0)
        return 0.5*x*(2.0*p.d()+x)/(p.d()*p.d())+0.5;
      if (x<p.d())
        return 0.5*x*(2.0*p.d()-x)/(p.d()*p.d())+0.5;
      return 1.0;
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<0.0 || x>1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      return icdf_(x);
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const tent_dist::param_type &p1, 
			 const tent_dist::param_type &p2) {
    return p1.m()==p2.m() && p1.d()==p2.d();
  }
  inline bool operator!=(const tent_dist::param_type &p1, 
			 const tent_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const tent_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.m() << ' ' << p.d()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     tent_dist::param_type &p) {
    double m, d;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> m >> utility::delim(' ')
       >> d >> utility::delim(')');
    if (in)
      p=tent_dist::param_type(m, d);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const tent_dist &g1, 
			 const tent_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const tent_dist &g1, 
			 const tent_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const tent_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[tent " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     tent_dist &g) {
    tent_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[tent ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
