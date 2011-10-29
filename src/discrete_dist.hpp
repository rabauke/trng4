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

#if !(defined TRNG_DISCRETE_DIST_HPP)

#define TRNG_DISCRETE_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <iomanip>
#include <istream>
#include <vector>

namespace trng {

  // non-uniform random number generator class
  class discrete_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      std::vector<double> P_;

      void calc_probabilities() {
        // build list with cumulative density function
        for (std::vector<double>::size_type i(1); i<P_.size(); ++i)
          P_[i]+=P_[i-1];
        for (std::vector<double>::size_type i(0); i<P_.size(); ++i)
          P_[i]/=P_.back();
      }

    public:
      template<typename iter>
      explicit param_type(iter first, iter last) :
	P_(first, last) {
	calc_probabilities();
      }
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
    explicit discrete_dist(iter first, iter last) : P(first, last) {
    }
    explicit discrete_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      return utility::discrete(utility::uniformco(r), 
			       P.P_.begin(), P.P_.end());
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      discrete_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return P.P_.size()-1; }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    // probability density function  
    double pdf(int x) const {
      if (x<0 || x>=P.P_.size())
        return 0.0;
      if (x==0)
        return P.P_[0];
      return P.P_[x]-P.P_[x-1];
    }
    // cumulative density function 
    double cdf(int x) const {
      if (x<0)
        return 0.0;
      if (x<P.P_.size())
        return P.P_[x];
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
    out << '(' << P.P_.size() << ' ';
    for (std::vector<double>::size_type i=0; i<P.P_.size(); ++i) {
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

