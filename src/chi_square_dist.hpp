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

#if !(defined TRNG_CHI_SQUARE_DIST_HPP)

#define TRNG_CHI_SQUARE_DIST_HPP

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
  class chi_square_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
    private:
      int nu_;
    public:
      int nu() const { return nu_; }
      void nu(int nu_new) { nu_=nu_new; }
      explicit param_type(int nu) : nu_(nu) {
      }
      friend class chi_square_dist;
    };
    
  private:
    param_type p;

    // inverse cumulative density function
    double icdf_(double x) const {
      if (x<=math::numeric_limits<double>::epsilon())
        return 0.0;
      const double kappa=0.5*p.nu();
      const double theta=2.0;
      if (kappa==1.0)  // special case of exponential distribution
        return -math::ln(1.0-x)*theta;
      const double ln_Gamma_kappa=math::ln_Gamma(kappa);
      double y=kappa, y_old;
      if (kappa<1.0 && x<0.5)
	y=x*x;
      int num_iterations=0;
      do {
	++num_iterations;
        y_old=y;
        double f0=math::GammaP(kappa, y)-x;
        double f1=math::pow(y, kappa-1.0)*math::exp(-y-ln_Gamma_kappa);
        double f2=f1*(kappa-1.0-y)/y;
        y-=f0/f1*(1+0.5*f0*f2/(f1*f1));
      } while (num_iterations<16 &&
               math::abs((y-y_old)/y)>16*math::numeric_limits<double>::epsilon());
      return y*theta;
    }

  public:
    // constructor
    explicit chi_square_dist(int nu) : p(nu) {
    }
    explicit chi_square_dist(const param_type &p) : p(p) {
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
      chi_square_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return 0; }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    int nu() const { return p.nu(); }
    void nu(int nu_new) { p.nu(nu_new); }
    // probability density function  
    double pdf(double x) const {
      if (x<0.0)
	return 0.0;
      else {
	x/=2.0;
	return math::pow(x, 0.5*p.nu()-1.0)/
	  (math::exp(x+math::ln_Gamma(0.5*p.nu()))*2.0);
      }
    }
    // cumulative density function 
    double cdf(double x) const {
      if (x<=0.0)
	return 0.0;
      else
	return math::GammaP(0.5*p.nu(), x/2.0);
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
  inline bool operator==(const chi_square_dist::param_type &p1, 
			 const chi_square_dist::param_type &p2) {
    return p1.nu()==p2.nu();
  }
  inline bool operator!=(const chi_square_dist::param_type &p1, 
			 const chi_square_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const chi_square_dist::param_type &p) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << p.nu()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     chi_square_dist::param_type &p) {
    int nu;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> nu >> utility::delim(')');
    if (in)
      p=chi_square_dist::param_type(nu);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const chi_square_dist &g1, 
			 const chi_square_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const chi_square_dist &g1, 
			 const chi_square_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const chi_square_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[chi_square " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     chi_square_dist &g) {
    chi_square_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[chi_square ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
