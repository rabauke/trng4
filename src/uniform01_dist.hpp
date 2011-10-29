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

#if !(defined TRNG_UNIFORM01_DIST_HPP)

#define TRNG_UNIFORM01_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  template<typename float_t=double>
  class uniform01_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    public:
      param_type() {
      }

      friend class uniform01_dist<float_t>;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
		 const param_type &p) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed |
		  std::ios_base::left);
	out << '('
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> utility::delim(')');
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    uniform01_dist() {
    }
    explicit uniform01_dist(const param_type &p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return utility::uniformco<result_type>(r);
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      return utility::uniformco<result_type>(r);
    }
    // property methods
    // min / max
    result_type min() const { return 0; }
    result_type max() const { return 1; }
    param_type param() const { return p; }
    void param(const param_type &p_new) { }
    // probability density function  
    result_type pdf(result_type x) const {
      if (x<0 or x>=1)
	return 0;
      return 1;
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      if (x<0)
	return 0;
      if (x>=1)
	return 1;
      return x;
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<0 or x>1) {
	errno=EDOM;
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      return x; 
    }
    
  };
  
  // -------------------------------------------------------------------

  // Equality comparable concept
  template<typename float_t>
  inline bool operator==(const typename uniform01_dist<float_t>::param_type &, 
			 const typename uniform01_dist<float_t>::param_type &) {
    return true;
  }

  template<typename float_t>
  inline bool operator!=(const typename uniform01_dist<float_t>::param_type &, 
			 const typename uniform01_dist<float_t>::param_type &) {
    return false;
  }
   
  // -------------------------------------------------------------------
  
  // Equality comparable concept
  template<typename float_t>
  inline bool operator==(const uniform01_dist<float_t> &g1, 
			 const uniform01_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const uniform01_dist<float_t> &g1, 
			 const uniform01_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform01_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[uniform01 " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform01_dist<float_t> &g) {
    typename uniform01_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[uniform01 ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
