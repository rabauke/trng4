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

#if !(defined TRNG_UNIFORM01_DIST_HPP)

#define TRNG_UNIFORM01_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class uniform01_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    public:
      explicit param_type() {
      }
      friend class uniform01_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit uniform01_dist() {
    }
    explicit uniform01_dist(const param_type &p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return utility::uniformco(r);
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      return utility::uniformco(r);
    }
    // property methods
    // min / max
    double min() const { return 0.0; }
    double max() const { return 1.0; }
    param_type param() const { return p; }
    void param(const param_type &p_new) { }
    // probability density function  
    double pdf(double x) const {
      if (x<0.0 || x>=1.0)
	return 0.0;
      return 1.0;
    }
    // cumulative density function 
    double cdf(double x) const {
      if (x<0.0)
	return 0;
      if (x>=1.0)
	return 1.0;
      return x;
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<0.0 || x>1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      return x; 
    }
    
    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &
    operator<<(std::basic_ostream<char_t, traits_t> &,
	       const uniform01_dist &);    
    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &
    operator>>(std::basic_istream<char_t, traits_t> &,
	       uniform01_dist &);
  };
  
  // -------------------------------------------------------------------

  // Equality comparable concept
  inline bool operator==(const uniform01_dist::param_type &, 
			 const uniform01_dist::param_type &) {
    return true;
  }
  inline bool operator!=(const uniform01_dist::param_type &, 
			 const uniform01_dist::param_type &) {
    return false;
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform01_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< ')';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform01_dist::param_type &p) {
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> utility::delim(')');
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------
  
  inline bool operator==(const uniform01_dist &g1, 
			 const uniform01_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const uniform01_dist &g1, 
			 const uniform01_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform01_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[uniform01 " << g.p << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform01_dist &g) {
    uniform01_dist::param_type p;
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
