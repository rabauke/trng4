// Copyright (C) 2006 Heiko Bauke <heiko.bauke@physik.uni-magdeburg.de>
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

#if !(defined TRNG_UNIFORM_DIST_HPP)

#define TRNG_UNIFORM_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class uniform_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double a_, b_, d_;
      double d() const { return d_; }
    public:
      double a() const { return a_; }
      void a(double a_new) { a_=a_new; d_=b_-a_; }
      double b() const { return b_; }
      void b(double b_new) { b_=b_new; d_=b_-a_; }
      explicit param_type(double a, double b) :
	a_(a), b_(b), d_(b-a) {
      }
      friend class uniform_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit uniform_dist(double a, double b) : p(a, b) {
    }
    explicit uniform_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return p.d()*utility::uniformco(r)+p.a();
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      uniform_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return p.a(); }
    double max() const { return p.b(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double a() const { return p.a(); }
    void a(double a_new) { p.a(a_new); }
    double b() const { return p.b(); }
    void b(double b_new) { p.b(b_new); }
    // probability density function  
    double pdf(double x) const {
      if (x<p.a() || x>=p.b())
	return 0.0;
      return 1.0/p.d();
    }
    // cumulative density function 
    double cdf(double x) const {
      if (x<p.a())
	return 0;
      if (x>=p.b())
	return 1.0;
      return (x-p.a())/p.d();
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<0.0 || x>1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      return x*p.d()+p.a(); 
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const uniform_dist::param_type &p1, 
			 const uniform_dist::param_type &p2) {
    return p1.a()==p2.a() && p1.b()==p2.b();
  }
  inline bool operator!=(const uniform_dist::param_type &p1, 
			 const uniform_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(16) << p.a() << ' ' << p.b()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform_dist::param_type &p) {
    double a, b;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> a >> utility::delim(' ')
       >> b >> utility::delim(')');
    if (in)
      p=uniform_dist::param_type(a, b);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const uniform_dist &g1, 
			 const uniform_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const uniform_dist &g1, 
			 const uniform_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[uniform " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform_dist &g) {
    uniform_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[uniform ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
