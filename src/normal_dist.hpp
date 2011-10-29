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

#if !(defined TRNG_NORMAL_DIST_HPP)

#define TRNG_NORMAL_DIST_HPP

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
  class normal_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      double mu_, sigma_;
    public:
      double mu() const { return mu_; }
      void mu(double mu_new) { mu_=mu_new; }
      double sigma() const { return sigma_; }
      void sigma(double sigma_new) { sigma_=sigma_new; }
      param_type(double mu, double sigma) : mu_(mu), sigma_(sigma) {
      }
      friend class normal_dist;
    };
    
  private:
    param_type p;
   
  public:
    // constructor
    normal_dist(double mu, double sigma) : p(mu, sigma) {
    }
    explicit normal_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      return icdf(utility::uniformoo(r));
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      normal_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return -math::numeric_limits<double>::infinity(); }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    double mu() const { return p.mu(); }
    void mu(double mu_new) { p.mu(mu_new); }
    double sigma() const { return p.sigma(); }
    void sigma(double sigma_new) { p.sigma(sigma_new); }
    // probability density function  
    double pdf(double x) const {
      double t(x-p.mu());
      return 0.39894228040143267794/p.sigma()*math::exp(-0.5*t*t/(p.sigma()*p.sigma()));
    }
    // cumulative density function 
    double cdf(double x) const {
      x-=p.mu();
      x/=p.sigma();
      return math::Phi(x);
    }
    // inverse cumulative density function 
    double icdf(double x) const {
      return math::inv_Phi(x)*p.sigma()+p.mu();
    }
  };
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const normal_dist::param_type &p1, 
			 const normal_dist::param_type &p2) {
    return p1.mu()==p2.mu() && p1.sigma()==p2.sigma();
  }
  inline bool operator!=(const normal_dist::param_type &p1, 
			 const normal_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const normal_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.mu() << ' ' << p.sigma() 
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     normal_dist::param_type &p) {
    double mu, sigma;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> mu >> utility::delim(' ')
       >> sigma >> utility::delim(')');
    if (in)
      p=normal_dist::param_type(mu, sigma);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const normal_dist &g1, 
			 const normal_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const normal_dist &g1, 
			 const normal_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const normal_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[normal " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     normal_dist &g) {
    normal_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[normal ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
