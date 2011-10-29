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

#if !(defined TRNG_WEIBULL_DIST_HPP)

#define TRNG_WEIBULL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  class weibull_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double theta_, beta_;
    public:
      double theta() const { return theta_; }
      void theta(double theta_new) { theta_=theta_new; }
      double beta() const { return beta_; }
      void beta(double beta_new) { beta_=beta_new; }
      explicit param_type(double theta, double beta) : theta_(theta), beta_(beta) {
      }
      friend class weibull_dist;
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    explicit weibull_dist(double theta, double beta) : p(theta, beta) {
    }
    explicit weibull_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return p.theta()*
	math::pow(-math::ln(utility::uniformoc(r)), 1.0/p.beta());
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      weibull_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0.0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double theta() const { return p.theta(); }
    void theta(double theta_new) { p.theta(theta_new); }
    double beta() const { return p.beta(); }
    void beta(double beta_new) { p.beta(beta_new); }
    // probability density function  
    double pdf(double x) const {
      if (x<0.0)
        return 0.0;
      x/=p.theta();
      if (x>0.0) {
        double t(math::pow(x, p.beta()));
        return p.beta()*t/x*math::exp(-t);
      }
      // x==0
      if (p.beta()==1.0)
        return 1.0/p.theta();
      if (p.beta()>1.0)
        return 0.0;
      return math::numeric_limits<double>::quiet_NaN();
    }
    // cumulative density function 
    double cdf(double x) const {
      x/=p.theta();
      if (x<=0.0)
        return 0.0;
      return 1.0-math::exp(-math::pow(x, p.beta()));
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      if (x<0.0 || x>=1.0) {
        errno=EDOM;
        return math::numeric_limits<double>::quiet_NaN();
      }
      return p.theta()*math::pow(-math::ln(1.0-x), 1.0/p.beta());
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const weibull_dist::param_type &p1, 
			 const weibull_dist::param_type &p2) {
    return p1.theta()==p2.theta() && p1.beta()==p2.beta();
  }
  inline bool operator!=(const weibull_dist::param_type &p1, 
			 const weibull_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const weibull_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.theta() << ' ' << p.beta() 
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     weibull_dist::param_type &p) {
    double theta, beta;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> theta >> utility::delim(' ')
       >> beta >> utility::delim(')');
    if (in)
      p=weibull_dist::param_type(theta, beta);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const weibull_dist &g1, 
			 const weibull_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const weibull_dist &g1, 
			 const weibull_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const weibull_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[weibull " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     weibull_dist &g) {
    weibull_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[weibull ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
