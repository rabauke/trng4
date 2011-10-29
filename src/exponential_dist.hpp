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

#if !(defined TRNG_EXPONENTIAL_DIST_HPP)

#define TRNG_EXPONENTIAL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class exponential_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double mu_;
    public:
      double mu() const { return mu_; }
      void mu(double mu_new) { mu_=mu_new; }
      explicit param_type(double mu) :	mu_(mu) {
      }
      friend class exponential_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit exponential_dist(double mu) : p(mu) {
    }
    explicit exponential_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return -p.mu()*math::ln(utility::uniformoc(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      exponential_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double mu() const { return p.mu(); }
    void mu(double mu_new) { p.mu(mu_new); }
    // probability density function  
    double pdf(double x) const {
      return x<0.0 ? 0.0 : math::exp(-x/p.mu())/p.mu();
    }
    // cumulative density function 
    double cdf(double x) const {
      return x<=0.0 ? 0.0 : 1.0-math::exp(-x/p.mu());
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<0.0 || x>1.0) {
	errno=EDOM;
	return math::numeric_limits<double>::quiet_NaN();
      }
      if (x==1.0)
        return math::numeric_limits<double>::infinity();
      return -p.mu()*math::ln(1.0-x);
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const exponential_dist::param_type &p1, 
			 const exponential_dist::param_type &p2) {
    return p1.mu()==p2.mu();
  }
  inline bool operator!=(const exponential_dist::param_type &p1, 
			 const exponential_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const exponential_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.mu() 
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     exponential_dist::param_type &p) {
    double mu;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> mu >> utility::delim(')');
    if (in)
      p=exponential_dist::param_type(mu);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const exponential_dist &g1, 
			 const exponential_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const exponential_dist &g1, 
			 const exponential_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const exponential_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[exponential " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     exponential_dist &g) {
    exponential_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[exponential ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
