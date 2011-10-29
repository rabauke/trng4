// Copyright (C) 2000-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#if !(defined TRNG_LOGNORMAL_DIST_HPP)

#define TRNG_LOGNORMAL_DIST_HPP

#include <trng/constants.hpp>
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
  class lognormal_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type mu_, sigma_;
    public:
      result_type mu() const { return mu_; }
      void mu(result_type mu_new) { mu_=mu_new; }
      result_type sigma() const { return sigma_; }
      void sigma(result_type sigma_new) { sigma_=sigma_new; }
      explicit param_type() : mu_(0), sigma_(1) {
      }
      param_type(result_type mu, result_type sigma) : mu_(mu), sigma_(sigma) {
      }

      friend class lognormal_dist;

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
            << p.mu() << ' ' << p.sigma() 
            << ')';
        out.flags(flags);
        return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
                 param_type &p) {
        float_t mu, sigma;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed |
                 std::ios_base::left);
        in >> utility::delim('(')
           >> mu >> utility::delim(' ')
           >> sigma >> utility::delim(')');
        if (in)
          p=param_type(mu, sigma);
        in.flags(flags);
        return in;
      }

    };
    
  private:
    param_type p;
    
  public:
    // constructor
    lognormal_dist(result_type mu, result_type sigma) : p(mu, sigma) {
    }
    explicit lognormal_dist(const param_type &p) : p(p) {
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
      lognormal_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return 0; }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type mu() const { return p.mu(); }
    void mu(result_type mu_new) { p.mu(mu_new); }
    result_type sigma() const { return p.sigma(); }
    void sigma(result_type sigma_new) { p.sigma(sigma_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      if (x<=0)
        return 0;
      result_type t((math::ln(x)-p.mu())/p.sigma());
      return math::constants<result_type>::one_over_sqrt_2pi()/(x*p.sigma())*math::exp(-t*t/2);
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      if (x<=0)
        return 0;
      return math::erfc(math::constants<result_type>::one_over_sqrt_2()*
			(p.mu()-math::ln(x))/p.sigma())/2;
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
    if (x<0 or x>1) {
        errno=EDOM;
        return math::numeric_limits<result_type>::quiet_NaN();
      }
      if (x==0)
        return 0;
      if (x==1)
        return math::numeric_limits<result_type>::infinity();
      return math::exp(math::inv_Phi(x)*p.sigma()+p.mu());
    }
  };
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename lognormal_dist<float_t>::param_type &p1, 
                         const typename lognormal_dist<float_t>::param_type &p2) {
    return p1.mu()==p2.mu() and p1.sigma()==p2.sigma();
  }

  template<typename float_t>
  inline bool operator!=(const typename lognormal_dist<float_t>::param_type &p1, 
                         const typename lognormal_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }  
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const lognormal_dist<float_t> &g1, 
                         const lognormal_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const lognormal_dist<float_t> &g1, 
                         const lognormal_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const lognormal_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[lognormal " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     lognormal_dist<float_t> &g) {
    typename lognormal_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[lognormal ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
