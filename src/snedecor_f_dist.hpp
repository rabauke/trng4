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
  template<typename float_t=double>
  class snedecor_f_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      int n_, m_;
    public:
      int n() const { return n_; }
      void n(int n_new) { n_=n_new; }
      int m() const { return m_; }
      void m(int m_new) { m_=m_new; }
      param_type() : n_(1), m_(1) {
      }
      param_type(int n, int m) : n_(n), m_(m) {
      }

      friend class snedecor_f_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
		 const param_type &p) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed |
		  std::ios_base::left);
	out << '('
	    << p.n() << ' ' << p.m()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
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

    };
    
  private:
    param_type p;

    // inverse cumulative density function
    result_type icdf_(result_type x) const {
      result_type t=math::inv_Beta_I(x, result_type(1)/result_type(2)*p.n(), result_type(1)/result_type(2)*p.m());
      return t/(1-t)*
	static_cast<result_type>(p.m())/static_cast<result_type>(p.n());
    }

  public:
    // constructor
    snedecor_f_dist(int n, int m) : p(n, m) {
    }
    explicit snedecor_f_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return icdf_(utility::uniformco<result_type>(r));
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      snedecor_f_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return 0; }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    int n() const { return p.n(); }
    void n(int n_new) { p.n(n_new); }
    int m() const { return p.m(); }
    void m(int m_new) { p.m(m_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      const result_type n=static_cast<result_type>(p.n()), m=static_cast<result_type>(p.m());
      // return math::exp(math::ln(n)*result_type(1)/result_type(2)*n + math::ln(m)*result_type(1)/result_type(2)*m - 
      // 		       math::ln(m+n*x)*result_type(1)/result_type(2)*(m+n) +
      // 		       math::ln(x)*(result_type(1)/result_type(2)*n-1))/math::Beta(result_type(1)/result_type(2)*n, result_type(1)/result_type(2)*m);
      return math::exp(result_type(1)/result_type(2)*math::ln(n/m)*n + 
		       math::ln(x)*(result_type(1)/result_type(2)*n-1) - 
		       math::ln(1+n*x/m)*(result_type(1)/result_type(2)*n + 
					    result_type(1)/result_type(2)*m) -
		       math::ln_Gamma(result_type(1)/result_type(2)*n) - 
		       math::ln_Gamma(result_type(1)/result_type(2)*m) + 
		       math::ln_Gamma(result_type(1)/result_type(2)*(n+m))
		       );
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      const result_type n=static_cast<result_type>(p.n()), 
	m=static_cast<result_type>(p.m());
      return math::Beta_I(n*x/(m+n*x), result_type(1)/result_type(2)*n, result_type(1)/result_type(2)*m);
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<=0 or x>=1) {
	errno=EDOM;
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      if (x==0)
	return 0;
      if (x==1)
	return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename snedecor_f_dist<float_t>::param_type &p1, 
			 const typename snedecor_f_dist<float_t>::param_type &p2) {
    return p1.n()==p2.n() and p1.m()==p2.m();
  }

  template<typename float_t>
  inline bool operator!=(const typename snedecor_f_dist<float_t>::param_type &p1, 
			 const typename snedecor_f_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------
  
  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const snedecor_f_dist<float_t> &g1, 
			 const snedecor_f_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }
  
  template<typename float_t>
  inline bool operator!=(const snedecor_f_dist<float_t> &g1, 
			 const snedecor_f_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const snedecor_f_dist<float_t> &g) {
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
	     snedecor_f_dist<float_t> &g) {
    typename snedecor_f_dist<float_t>::param_type p;
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
