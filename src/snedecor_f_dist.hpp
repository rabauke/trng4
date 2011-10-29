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

#if !(defined TRNG_SNEDECOR_F_DIST_HPP)

#define TRNG_SNEDECOR_F_DIST_HPP

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
  class snedecor_f_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      int n_, m_;
    public:
      int n() const { return n_; }
      void n(int n_new) { n_=n_new; }
      int m() const { return m_; }
      void m(int m_new) { m_=m_new; }
      explicit param_type(int n, int m) : n_(n), m_(m) {
      }
      friend class snedecor_f_dist;
    };
    
  private:
    param_type p;

    // inverse cumulative density function
    double icdf_(double x) const {
      double t=math::inv_Beta_I(x, 0.5*p.n(), 0.5*p.m());
      return t/(1.0-t)*
	static_cast<double>(p.m())/static_cast<double>(p.n());
    }

  public:
    // constructor
    explicit snedecor_f_dist(int n, int m) : p(n, m) {
    }
    explicit snedecor_f_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return icdf_(utility::uniformco(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      snedecor_f_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0.0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    int n() const { return p.n(); }
    void n(int n_new) { p.n(n_new); }
    int m() const { return p.m(); }
    void m(int m_new) { p.m(m_new); }
    // probability density function  
    double pdf(double x) const {
      const double n=static_cast<double>(p.n()), m=static_cast<double>(p.m());
      // return math::exp(math::ln(n)*0.5*n + math::ln(m)*0.5*m - 
      // 		       math::ln(m+n*x)*0.5*(m+n) +
      // 		       math::ln(x)*(0.5*n-1.0))/math::Beta(0.5*n, 0.5*m);
      return math::exp(0.5*math::ln(n/m)*n + 
		       math::ln(x)*(0.5*n-1.0) - math::ln(1.0+n*x/m)*(0.5*n+0.5*m)
		       -math::ln_Gamma(0.5*n)-math::ln_Gamma(0.5*m)+math::ln_Gamma(0.5*(n+m))
		       );
    }
    // cumulative density function 
    double cdf(double x) const {
      const double n=static_cast<double>(p.n()), m=static_cast<double>(p.m());
      return math::Beta_I(n*x/(m+n*x), 0.5*n, 0.5*m);
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
      return icdf_(x);
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const snedecor_f_dist::param_type &p1, 
			 const snedecor_f_dist::param_type &p2) {
    return p1.n()==p2.n() && p1.m()==p2.m();
  }
  inline bool operator!=(const snedecor_f_dist::param_type &p1, 
			 const snedecor_f_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const snedecor_f_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.n() << ' ' << std::setprecision(17) << p.m()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     snedecor_f_dist::param_type &p) {
    int n, m;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> n >> utility::delim(' ')
       >> m >> utility::delim(')');
    if (in)
      p=snedecor_f_dist::param_type(n, m);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const snedecor_f_dist &g1, 
			 const snedecor_f_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const snedecor_f_dist &g1, 
			 const snedecor_f_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const snedecor_f_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[snedecor_f " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     snedecor_f_dist &g) {
    snedecor_f_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[snedecor_f ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
