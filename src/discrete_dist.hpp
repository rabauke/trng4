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

#if !(defined TRNG_DISCRETE_DIST_HPP)

#define TRNG_DISCRETE_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <iomanip>
#include <istream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace trng {

  namespace detail {

    template<typename T>
    inline T log2_floor(T x) {
      T y(0);
      while (x>0) {
	x>>=1;
	++y;
      };
      --y;
      return y;
    }
    
    template<typename T>
    inline T log2_ceil(T x) {
      T y(log2_floor(x));
      if ((T(1)<<y)<x)
        ++y;
      return y;
    }
    
    template<typename T>
    inline T pow2(T x) {
      return T(1) << x;
    }
    
  }

  // non-uniform random number generator class
  class discrete_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      typedef std::vector<double>::size_type size_type;
      std::vector<double> P_;
      size_type N, offset, layers;
      
    public:
      template<typename iter>
      param_type(iter first, iter last) :
	P_(first, last) {
	N=P_.size();
	layers=detail::log2_ceil(N);
	offset=detail::pow2(layers)-1;
	P_.resize(N+offset);
	std::copy_backward(P_.begin(), P_.begin()+N, P_.end());
	std::fill(P_.begin(), P_.begin()+N, 0);
	update_all_layers();
      }
      explicit param_type(int n) :
	P_(n, 1.0) {
	N=P_.size();
	layers=detail::log2_ceil(N);
	offset=detail::pow2(layers)-1;
	P_.resize(N+offset);
	std::copy_backward(P_.begin(), P_.begin()+N, P_.end());
	std::fill(P_.begin(), P_.begin()+N, 0);
	update_all_layers();
      }
    private:
      void update_layer(size_type layer, size_type n) {
	size_type first=detail::pow2(layer)-1, last=first+n;
	for (size_type i=first; i<last; ++i, ++i) 
	  if (i+1<last)
	    P_[(i-1)/2]=P_[i]+P_[i+1];
	  else
	    P_[(i-1)/2]=P_[i];
      }
      void update_all_layers() {
	size_type layer=layers;
	if (layer>0) {
	  update_layer(layer, N);
	  --layer;
	}
	while (layer>0) {
	  update_layer(layer, detail::pow2(layer));
	  --layer;
	}
      }
    public:
      friend class discrete_dist;
      friend bool operator==(const param_type &, const param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &,
		 const discrete_dist::param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &,
		 discrete_dist::param_type &);
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    template<typename iter>
    discrete_dist(iter first, iter last) : P(first, last) {
    }
    explicit discrete_dist(int N) : P(N) {
    }
    explicit discrete_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      double u=utility::uniformco(r)*P.P_[0];
      int x=0;
      while (x<P.offset) {
	if (u<P.P_[2*x+1]) {
	  x=2*x+1;
	} else {
	  u-=P.P_[2*x+1];
	  x=2*x+2;
	}
      }
      return x-P.offset;
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      discrete_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return P.N-1; }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    void param(int x, double p) {
      x+=P.offset;
      P.P_[x]=p;
      if (x>0) {
	do {
	  x=(x-1)/2;
	  P.P_[x]=P.P_[2*x+1]+P.P_[2*x+2];
	} while (x>0);
      }
    }
    // probability density function  
    double pdf(int x) const {
      return (x<0 || x>=static_cast<int>(P.N)) ? 0.0 : P.P_[x+P.offset]/P.P_[0];
    }
    // cumulative density function 
    double cdf(int x) const {
      if (x<0)
	return 0.0;
      if (x<static_cast<int>(P.N))
	return std::accumulate(&P.P_[P.offset], &P.P_[x+P.offset+1], 0.0)/P.P_[0];
      return 1.0;
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const discrete_dist::param_type &p1, 
			 const discrete_dist::param_type &p2) {
    return p1.P_==p2.P_;
  }
  inline bool operator!=(const discrete_dist::param_type &p1, 
			 const discrete_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const discrete_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '(' << P.N << ' ';
    for (std::vector<double>::size_type i=P.offset; i<P.P_.size(); ++i) {
      out << std::setprecision(17) << P.P_[i];
      if (i+1<P.P_.size())
	out << ' ';
    }
    out	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     discrete_dist::param_type &P) {
    double p;
    std::vector<double>::size_type n;
    std::vector<double> P_new;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> n >> utility::delim(' ');
    for (std::vector<double>::size_type i=0; i<n; ++i) {
      in >> p;
      if (i+1<n)
	in >> utility::delim(' ');
      P_new.push_back(p);
    }
    in >> utility::delim(')');
    if (in)
      P=discrete_dist::param_type(P_new.begin(), P_new.end());
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const discrete_dist &g1, 
			 const discrete_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const discrete_dist &g1, 
			 const discrete_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const discrete_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[discrete " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     discrete_dist &g) {
    discrete_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[discrete ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
