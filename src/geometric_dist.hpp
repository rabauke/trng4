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

#if !(defined TRNG_GEOMETRIC_DIST_HPP)

#define TRNG_GEOMETRIC_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>

namespace trng {

  // non-uniform random number generator class
  class geometric_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      double p_, q_, one_over_ln_q_;
     
      double q() const {  return q_;  }
      double one_over_ln_q() const {  return one_over_ln_q_;  }

    public:
      double p() const { return p_; }
      void p(double p_new) { 
	p_=p_new;  q_=1.0-p_;  one_over_ln_q_=1.0/math::ln(q_);
      }
      explicit param_type(double p) :
	p_(p), q_(1.0-p_), one_over_ln_q_(1.0/math::ln(q_)) {
      }
      friend class geometric_dist;
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    explicit geometric_dist(double p) : P(p) {
    }
    explicit geometric_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      return static_cast<int>(math::ln(utility::uniformoo(r))*
			      P.one_over_ln_q());
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      geometric_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return math::numeric_limits<int>::max(); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    double p() const { return P.p(); }
    void p(double p_new) { P.p(p_new); }
    // probability density function  
    double pdf(int x) const {
      return x<0 ? 0.0 : P.p()*math::pow(P.q(), static_cast<double>(x));
    }
    // cumulative density function 
    double cdf(int x) const {
      return x<0 ? 0.0 : 1.0-math::pow(P.q(), static_cast<double>(x+1));
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const geometric_dist::param_type &p1, 
			 const geometric_dist::param_type &p2) {
    return p1.p()==p2.p();
  }
  inline bool operator!=(const geometric_dist::param_type &p1, 
			 const geometric_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const geometric_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(16) << P.p()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     geometric_dist::param_type &P) {
    double p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> p >> utility::delim(')');
    if (in)
      P=geometric_dist::param_type(p);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const geometric_dist &g1, 
			 const geometric_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const geometric_dist &g1, 
			 const geometric_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const geometric_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[geometric " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     geometric_dist &g) {
    geometric_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[geometric ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif

