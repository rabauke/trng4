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

#if !(defined TRNG_UNIFORM_DIST_HPP)

#define TRNG_UNIFORM_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>

namespace trng {

  // uniform random number generator class
  template<typename float_t=double>
  class uniform_dist {
  public:
    typedef float_t result_type;
    class param_type;
    
    class param_type {
    private:
      result_type a_, b_, d_;
      result_type d() const { return d_; }
    public:
      result_type a() const { return a_; }
      void a(result_type a_new) { a_=a_new; d_=b_-a_; }
      result_type b() const { return b_; }
      void b(result_type b_new) { b_=b_new; d_=b_-a_; }
      param_type() :
	a_(0), b_(1), d_(1) {
      }
      param_type(result_type a, result_type b) :
	a_(a), b_(b), d_(b-a) {
      }

      friend class uniform_dist<float_t>;

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
	    << p.a() << ' ' << p.b()
	    << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 param_type &p) {
	float_t a, b;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> a >> utility::delim(' ')
	   >> b >> utility::delim(')');
	if (in)
	  p=param_type(a, b);
	in.flags(flags);
	return in;
      }
      
    };
    
  private:
    param_type p;
    
  public:
    // constructor
    uniform_dist(result_type a, result_type b) : p(a, b) {
    }
    explicit uniform_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      return p.d()*utility::uniformco<result_type>(r)+p.a();
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      uniform_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return p.a(); }
    result_type max() const { return p.b(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
    result_type a() const { return p.a(); }
    void a(result_type a_new) { p.a(a_new); }
    result_type b() const { return p.b(); }
    void b(result_type b_new) { p.b(b_new); }
    // probability density function  
    result_type pdf(result_type x) const {
      if (x<p.a() or x>=p.b())
	return 0.0;
      return 1.0/p.d();
    }
    // cumulative density function 
    result_type cdf(result_type x) const {
      if (x<p.a())
	return 0;
      if (x>=p.b())
	return 1;
      return (x-p.a())/p.d();
    }
    // inverse cumulative density function 
    result_type icdf(result_type x) const {
      if (x<0 or x>1) {
	errno=EDOM;
	return math::numeric_limits<result_type>::quiet_NaN();
      }
      return x*p.d()+p.a(); 
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const typename uniform_dist<float_t>::param_type &p1, 
			 const typename uniform_dist<float_t>::param_type &p2) {
    return p1.a()==p2.a() and p1.b()==p2.b();
  }

  template<typename float_t>
  inline bool operator!=(const typename uniform_dist<float_t>::param_type &p1, 
			 const typename uniform_dist<float_t>::param_type &p2) {
    return not (p1==p2);
  }
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const uniform_dist<float_t> &g1, 
			 const uniform_dist<float_t> &g2) {
    return g1.param()==g2.param();
  }

  template<typename float_t>
  inline bool operator!=(const uniform_dist<float_t> &g1, 
			 const uniform_dist<float_t> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const uniform_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[uniform " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     uniform_dist<float_t> &g) {
    typename uniform_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[uniform ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
