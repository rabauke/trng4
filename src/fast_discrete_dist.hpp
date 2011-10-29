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

#if !(defined TRNG_FAST_DISCRETE_DIST_HPP)

#define TRNG_FAST_DISCRETE_DIST_HPP

// Algorithm described in
//
// Richard A. Kronmal; Arthur V. Peterson, Jr.
// On the Alias Method for Generating Random Variables from a Discrete Distribution
// The American Statistician, Vol. 33, No. 4. (Nov., 1979), pp. 214-218.
// 
// http://links.jstor.org/sici?sici=0003-1305%28197911%2933%3A4%3C214%3AOTAMFG%3E2.0.CO%3B2-1

#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <iomanip>
#include <istream>
#include <vector>
#include <numeric>
#include <functional>

namespace trng {

  // non-uniform random number generator class
  class fast_discrete_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      typedef std::vector<double>::size_type size_type;
      std::vector<double> P, F;
      std::vector<int> L;
      size_type N;
      
    public:
      template<typename iter>
      param_type(iter first, iter last) :
	P(first, last), F(P.size()), L(P.size()), N(P.size()) {
	update();
      }
      explicit param_type(int n) :
	P(n, 1.0), F(P.size()), L(P.size()), N(P.size()) {
	update();
      }
    private:
      void update() {
	double s=std::accumulate(P.begin(), P.end(), 0.0);
	if (s>0.0) {
	  for (int i=0; i<P.size(); ++i) 
	    P[i]/=s;
	  std::vector<int> G, S;
	  G.reserve(N);
	  S.reserve(N);
	  for (int i=0; i<F.size(); ++i) {
	    F[i]=N*P[i];
	    if (F[i]<1.0)
	      S.push_back(i);
	    else
	      G.push_back(i);
	  }
	  while (not S.empty()) {
	    int k=G.back(), j=S.back();
	    L[j]=k;
	    F[k]-=1.0-F[j];
	    S.pop_back();
	    if (F[k]<1.0) {
	      G.pop_back();
	      S.push_back(k);
	    }
	  }
	}
      }
    public:
      friend class fast_discrete_dist;
      friend bool operator==(const param_type &, const param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &
      operator<<(std::basic_ostream<char_t, traits_t> &,
		 const fast_discrete_dist::param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &
      operator>>(std::basic_istream<char_t, traits_t> &,
		 fast_discrete_dist::param_type &);
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    template<typename iter>
    fast_discrete_dist(iter first, iter last) : P(first, last) {
    }
    explicit fast_discrete_dist(int N) : P(N) {
    }
    explicit fast_discrete_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      double U=utility::uniformco<double>(r)*P.N;
      int I=static_cast<int>(U);
      return U-I<=P.F[I] ? I : P.L[I];
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      fast_discrete_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return P.N-1; }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    // probability density function  
    double pdf(int x) const {
      return (x<0 or x>=static_cast<int>(P.N)) ? 0.0 : P.P[x];
    }
    // cumulative density function 
    double cdf(int x) const {
      if (x<0)
	return 0.0;
      if (x<static_cast<int>(P.N))
	return std::accumulate(P.P.begin(), P.P.begin()+x+1, 0.0);
      return 1.0;
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const fast_discrete_dist::param_type &p1, 
			 const fast_discrete_dist::param_type &p2) {
    return p1.P==p2.P;
  }
  inline bool operator!=(const fast_discrete_dist::param_type &p1, 
			 const fast_discrete_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const fast_discrete_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '(' << P.N << ' ';
    for (std::vector<double>::size_type i=0; i<P.P.size(); ++i) {
      out << std::setprecision(17) << P.P[i];
      if (i+1<P.P.size())
	out << ' ';
    }
    out	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     fast_discrete_dist::param_type &P) {
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
      P=fast_discrete_dist::param_type(P_new.begin(), P_new.end());
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const fast_discrete_dist &g1, 
			 const fast_discrete_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const fast_discrete_dist &g1, 
			 const fast_discrete_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const fast_discrete_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[fast_discrete " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     fast_discrete_dist &g) {
    fast_discrete_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[fast_discrete ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif
