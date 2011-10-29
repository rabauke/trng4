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

#if !(defined TRNG_STUDENT_T_DIST_HPP)

#define TRNG_STUDENT_T_DIST_HPP

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
  class student_t_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      int nu_;
    public:
      int nu() const { return nu_; }
      void nu(int nu_new) { nu_=nu_new; }
      param_type() : nu_(1) {
      }
      explicit param_type(int nu) : nu_(nu) {
      }

      friend class student_t_dist;
      
      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
                 const param_type &p) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed |
                  std::ios_base::left);
        out << '('
            << p.nu()
            << ')';
        out.flags(flags);
        return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
                 param_type &p) {
        int nu;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed |
                 std::ios_base::left);
        in >> utility::delim('(')
           >> nu >> utility::delim(')');
        if (in)
          p=param_type(nu);
        in.flags(flags);
        return in;
      }

  };
    
  private:
    param_type p;

    // inverse cumulative density function
    result_type icdf_(result_type x) const {
      result_type t=math::inv_Beta_I(x, p.nu()/result_type(2), p.nu()/result_type(2));
      return math::sqrt(p.nu()/(t*(1-t)))*(t-result_type(1)/result_type(2));
    }

  public:
    // constructor
    explicit student_t_dist(int nu) : p(nu) {
    }
    explicit student_t_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return icdf_(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      student_t_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return -math::numeric_limits<result_type>::infinity(); }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    int nu() const { return p.nu(); }
    void nu(int nu_new) { p.nu(nu_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      result_type norm=
	math::exp(math::ln_Gamma((p.nu()+1)/result_type(2))-
		  math::ln_Gamma(p.nu()/result_type(2)))/
	math::sqrt(math::constants<result_type>::pi()*p.nu());
      return norm*math::pow(1+x*x/p.nu(), (p.nu()+1)/result_type(-2));
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      result_type t1=+math::sqrt(x*x+p.nu());
      result_type t2=(x+t1)/(2*t1);
      return math::Beta_I(t2, p.nu()/result_type(2), p.nu()/result_type(2));
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<=0 or x>=1) {
	errno=EDOM;
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      if (x==0)
	return -math::numeric_limits<result_type>::infinity();
      if (x==1)
	return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename student_t_dist<float_t>::param_type &p1, 
			 const typename student_t_dist<float_t>::param_type &p2) {
    return p1.nu()==p2.nu();
  }

  template<typename float_t>
  inline bool operator!=(const typename student_t_dist<float_t>::param_type &p1, 
			 const typename student_t_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const student_t_dist<float_t> &g1, 
			 const student_t_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const student_t_dist<float_t> &g1, 
			 const student_t_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const student_t_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[student_t " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     student_t_dist<float_t> &g) {
    typename student_t_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[student_t ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
