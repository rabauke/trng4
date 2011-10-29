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

#if !(defined TRNG_BERNOULLI_DIST_HPP)

#define TRNG_BERNOULLI_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>

namespace trng {

  // non-uniform random number generator class
  template<typename T>
  class bernoulli_dist {
  public:
    typedef T result_type;
    class param_type;
    
    class param_type {
    private:
      double p_;
      T head_, tail_;
    public:
      double p() const { return p_; }
      void p(double p_new) { p_=p_new; }
      T head() const { return head_; }
      void head(const T &head_new) { head_=head_new; }
      T tail() const { return tail_; }
      void tail(const T &tail_new) { tail_=tail_new; }
      param_type(double p, const T &head, const T &tail) :
	p_(p), head_(head), tail_(tail) {
      }
      friend class bernoulli_dist;
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    bernoulli_dist(double p, const T &head, const T &tail) : 
      P(p, head, tail) {
    }
    explicit bernoulli_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    T operator()(R &r) {
      return utility::uniformco<double>(r)<P.p() ? P.head() : P.tail();
    }
    template<typename R>
    T operator()(R &r, const param_type &P) {
      bernoulli_dist g(P);
      return g(r);
    }
    // property methods
    T min() const { return P.head(); }
    T max() const { return P.tail(); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    double p() const { return P.p(); }
    void p(double p_new) { P.p(p_new); }
    T head() const { return P.head(); }
    void head(const T & head_new) { P.head(head_new); }
    T tail() const { return P.tail(); }
    void tail(const T & tail_new) { P.tail(tail_new); }
    // probability density function  
    double pdf(const T &x) const {
      if (x==P.head())
	return P.p();
      else if (x==P.tail())
	return 1.0-P.p();
      return 0.0;
    }
    // cumulative density function 
    double cdf(const T &x) const {
      if (x==P.head())
	return P.p();
      else if (x==P.tail())
	return 1.0;
      return 0.0;
    }
    
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename T>
  inline bool operator==(const typename bernoulli_dist<T>::param_type &p1, 
			 const typename bernoulli_dist<T>::param_type &p2) {
    return p1.p()==p2.p() and p1.head()==p2.head() and p1.tail()==p2.tail();
  }
  template<typename T>
  inline bool operator!=(const typename bernoulli_dist<T>::param_type &p1, 
			 const typename bernoulli_dist<T>::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename T>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const typename bernoulli_dist<T>::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< P.p() << ' ' << P.head() << ' ' << P.tail()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename T>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     typename bernoulli_dist<T>::param_type &P) {
    double p;
    T head, tail;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> p >> utility::delim(' ')
       >> head >> utility::delim(' ')
       >> tail >> utility::delim(')');
    if (in)
      P=typename bernoulli_dist<T>::param_type(p, head, tail);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------
  
  // EqualityComparable concept
  template<typename T>
  inline bool operator==(const bernoulli_dist<T> &g1, 
			 const bernoulli_dist<T> &g2) {
    return g1.param()==g2.param();
  }
  template<typename T>
  inline bool operator!=(const bernoulli_dist<T> &g1, 
			 const bernoulli_dist<T> &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t, typename T>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const bernoulli_dist<T> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[bernoulli " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t, typename T>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     bernoulli_dist<T> &g) {
    typename bernoulli_dist<T>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[bernoulli ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif

