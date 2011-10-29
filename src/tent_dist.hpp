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
  template<typename float_t=double>
  class tent_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type m_, d_;
    public:
      result_type m() const { return m_; }
      void m(result_type m_new) { m_=m_new; }
      result_type d() const { return d_; }
      void d(result_type d_new) { d_=d_new; }
      param_type() : m_(0), d_(1) {
      }
      param_type(result_type m, result_type d) : m_(m), d_(d) {
      }

      friend class tent_dist;

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
	    << p.m() << ' ' << p.d()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t m, d;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> m >> utility::delim(' ')
	   >> d >> utility::delim(')');
	if (in)
	  p=param_type(m, d);
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
    result_type icdf_(result_type x) const {
      return x<result_type(1)/result_type(2) ?
        (-1+math::sqrt(    2*x))*p.d()+p.m() :
        ( 1-math::sqrt(2-2*x))*p.d()+p.m();
    }

  public:
    // constructor
    tent_dist(result_type m, result_type d) : p(m, d) {
    }
    explicit tent_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return icdf_(utility::uniformcc<result_type>(r));
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      tent_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return p.m()-p.d(); }
    result_type max() const { return p.m()+p.d(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type m() const { return p.m(); }
    void m(result_type m_new) { p.m(m_new); }
    result_type d() const { return p.d(); }
    void d(result_type d_new) { p.d(d_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      x-=p.m();
      if (x<=-p.d() or x>=p.d())
        return 0;
      return x<0 ? (1+x/p.d())/p.d() : (1-x/p.d())/p.d();
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      x-=p.m();
      if (x<=-p.d())
        return 0;
      if (x<=0)
        return (x+p.d())*(x+p.d())/(result_type(2)*p.d()*p.d());
      if (x<p.d())
        return 1-(x-p.d())*(x-p.d())/(result_type(2)*p.d()*p.d());
      return 1;
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<0 or x>1) {
	errno=EDOM;
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      return icdf_(x);
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename tent_dist<float_t>::param_type &p1, 
			 const typename tent_dist<float_t>::param_type &p2) {
    return p1.m()==p2.m() and p1.d()==p2.d();
  }

  template<typename float_t>
  inline bool operator!=(const typename tent_dist<float_t>::param_type &p1, 
			 const typename tent_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const tent_dist<float_t> &g1, 
			 const tent_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const tent_dist<float_t> &g1, 
			 const tent_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const tent_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[tent " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     tent_dist<float_t> &g) {
    typename tent_dist<float_t>::param_type p;
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
