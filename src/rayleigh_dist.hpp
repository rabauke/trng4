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
  template<typename float_t=double>
  class rayleigh_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      result_type nu_;
    public:
      result_type nu() const { return nu_; }
      void nu(result_type nu_new) { nu_=nu_new; }
      param_type() : nu_(1) {
      }
      explicit param_type(result_type nu) : nu_(nu) {
      }

      friend class rayleigh_dist;

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
	    << p.nu()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t nu;
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

    };
    
  private:
    param_type p;
   
  public:
    // constructor
    explicit rayleigh_dist(result_type nu) : p(nu) {
    }
    explicit rayleigh_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return icdf(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      rayleigh_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return 0; }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type nu() const { return p.nu(); }
    void nu(result_type nu_new) { p.nu(nu_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      if (x<=0)
	return 0;
      result_type t=x/(p.nu()*p.nu());
      return t*math::exp(-t*x/2);
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      if (x<=0)
	return 0;
      return 1-math::exp(-x*x/(2*p.nu()*p.nu()));
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      return p.nu()*math::sqrt(-2*math::ln(1-x));
    }
  };
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename rayleigh_dist<float_t>::param_type &p1, 
			 const typename rayleigh_dist<float_t>::param_type &p2) {
    return p1.nu()==p2.nu();
  }

  template<typename float_t>
  inline bool operator!=(const typename rayleigh_dist<float_t>::param_type &p1, 
			 const typename rayleigh_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const rayleigh_dist<float_t> &g1, 
			 const rayleigh_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const rayleigh_dist<float_t> &g1, 
			 const rayleigh_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const rayleigh_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[rayleigh " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     rayleigh_dist<float_t> &g) {
    typename rayleigh_dist<float_t>::param_type p;
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
