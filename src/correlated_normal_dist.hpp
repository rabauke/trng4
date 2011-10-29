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

#if !(defined TRNG_CORRELATED_NORMAL_DIST_HPP)

#define TRNG_CORRELATED_NORMAL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>
#include <vector>
#include <iterator>
#include <algorithm>

namespace trng {

  // uniform random number generator class
  class correlated_normal_dist {
  public:
    typedef double result_type;
    class param_type;
    
    class param_type {
      typedef std::vector<double>::size_type size_type;
      std::vector<double> H_;
      size_type d_;

      void Cholesky_factorization() {
	for (size_type i=0; i<d_; ++i) {
	  double t;
	  for (size_type k=0; k<i; ++k) {
	    t=0.0;
	    for (size_type j=0; j<k; ++j)
	      t+=H_[i*d_ + j]*H_[k*d_ + j];
	    H_[i*d_ + k]=( H_[i*d_ + k] - t )/H_[k*d_ + k];
	  }
	  t=0.0;
	  for (size_type j=0; j<i; ++j)
	    t+=H_[i*d_ + j]*H_[i*d_ + j];
	  H_[i*d_ + i]=trng::math::sqrt(H_[i*d_ + i]-t);
	  for (size_type k=i+1; k<d_; ++k)
	    H_[i*d_ + k]=0;
	}
      }
      
      double H_times(const std::vector<double> &normal) {
	size_type d(normal.size());
	double y=0.0;
	for (size_type j=0; j<d; ++j)
	  y+=H_[(d-1)*d_ + j]*normal[j];
	return y;
      }

      size_type d() const {
	return d_;
      }

    public:
      template<typename iter>
      param_type(iter first, iter last) {
	d_=static_cast<size_type>(trng::math::sqrt(static_cast<double>(std::distance(first, last))));
	H_.reserve(d_*d_);
	for (size_type i=0; i<d_*d_; ++i, ++first)
	  H_.push_back(*first);
	Cholesky_factorization();
      }

      friend class correlated_normal_dist;

      // EqualityComparable concept
      friend inline bool operator==(const correlated_normal_dist::param_type &p1, 
				    const correlated_normal_dist::param_type &p2) {
	return p1.H_==p2.H_;
      }
      friend inline bool operator!=(const correlated_normal_dist::param_type &p1, 
				    const correlated_normal_dist::param_type &p2) {
	return !(p1==p2);
      }

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &out,
		 const correlated_normal_dist::param_type &p) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed |
		  std::ios_base::left);
	out << '(' << p.d_ << std::setprecision(17); 
	for (unsigned int i=0; i<p.d_*p.d_; ++i) 
	  out << ' ' << p.H_[i];
	out << ')';
	out.flags(flags);
	return out;
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &in,
		 correlated_normal_dist::param_type &p) {
	std::vector<double> H, normal;
	unsigned int d;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed |
		 std::ios_base::left);
 	in >> utility::delim('(') >> d ;
	for (unsigned int i=0; i<d*d; ++i) {
	  double t=0;
	  in >> utility::delim(' ') >> t;
	  H.push_back(t);
	}
	in >> utility::delim(')');
 	if (in) {
	  p.d_=d;
	  p.H_=H;
	}
	in.flags(flags);
	return in;
      }
    };
    
  private:
    param_type p;
    std::vector<double> normal;

  public:
    // constructor
    template<typename iter>
    correlated_normal_dist(iter first, iter last) : p(first, last) {
    }
    explicit correlated_normal_dist(const param_type &p) : p(p) {
    }
    // reset internal state
    void reset() { normal.clear(); }
    // random numbers
    template<typename R>
    double operator()(R &r) {
      normal.push_back(trng::math::inv_Phi(utility::uniformoo(r)));
      double y=p.H_times(normal);
      if (normal.size()==p.d())
	normal.clear();
      return y;
    }
    template<typename R>
    double operator()(R &r, const param_type &p) {
      correlated_normal_dist g(p);
      return g(r);
    }
    // property methods
    double min() const { return -math::numeric_limits<double>::infinity(); }
    double max() const { return math::numeric_limits<double>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p=p_new; }
  };
    
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const correlated_normal_dist &g1, 
			 const correlated_normal_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const correlated_normal_dist &g1, 
			 const correlated_normal_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const correlated_normal_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[correlated_normal " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     correlated_normal_dist &g) {
    correlated_normal_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[correlated_normal ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
