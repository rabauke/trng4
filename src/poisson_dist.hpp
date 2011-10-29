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

#if !(defined TRNG_POISSON_DIST_HPP)

#define TRNG_POISSON_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>

namespace trng {

  // non-uniform random number generator class
  class poisson_dist {
  public:
    typedef int result_type;
    class param_type;
    
    class param_type {
    private:
      double mu_;
      std::vector<double> P_;
      
      void calc_probabilities() {
	P_=std::vector<double>();
	int x=0;
	double p=0.0;
	while (p<1.0-1.0/4096.0) {
	  p=math::GammaQ(x+1.0, mu_);
	  P_.push_back(p);
	  ++x;
	}
      }
      
    public:
      double mu() const { return mu_; }
      void mu(double mu_new) { mu_=mu_new;  calc_probabilities(); }
      explicit param_type(double mu) : mu_(mu) {
	calc_probabilities();
      }
      friend class poisson_dist;
    };
    
  private:
    param_type P;
    
  public:
    // constructor
    explicit poisson_dist(double mu) : P(mu) {
    }
    explicit poisson_dist(const param_type &P) : P(P) {
    }
    // reset internal state
    void reset() { }
    // random numbers
    template<typename R>
    int operator()(R &r) {
      double p(utility::uniformco<double>(r));
      int x(utility::discrete(p, P.P_.begin(), P.P_.end()));
      if (x+1==P.P_.size())
        while (p>cdf(x))
          ++x;
      return x;
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      poisson_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return math::numeric_limits<int>::max(); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P=P_new; }
    double mu() const { return P.mu(); }
    void mu(double mu_new) { P.mu(mu_new); }
    // probability density function  
    double pdf(int x) const {
      return x<0 ? 0.0 : math::exp(-P.mu()
				   -math::ln_Gamma(x+1.0)
				   +x*math::ln(P.mu()));
    }
    // cumulative density function 
    double cdf(int x) const {
      return x<0 ? 0.0 : math::GammaQ(x+1.0, P.mu());
    }
  };
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const poisson_dist::param_type &p1, 
			 const poisson_dist::param_type &p2) {
    return p1.mu()==p2.mu();
  }
  inline bool operator!=(const poisson_dist::param_type &p1, 
			 const poisson_dist::param_type &p2) {
    return !(p1==p2);
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const poisson_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << '('
	<< std::setprecision(17) << P.mu()
	<< ')';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     poisson_dist::param_type &P) {
    double mu;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::delim('(')
       >> mu >> utility::delim(')');
    if (in)
      P=poisson_dist::param_type(mu);
    in.flags(flags);
    return in;
  }
  
  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const poisson_dist &g1, 
			 const poisson_dist &g2) {
    return g1.param()==g2.param();
  }
  inline bool operator!=(const poisson_dist &g1, 
			 const poisson_dist &g2) {
    return g1.param()!=g2.param();
  }
  
  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &
  operator<<(std::basic_ostream<char_t, traits_t> &out,
	     const poisson_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed |
	      std::ios_base::left);
    out << "[poisson " << g.param() << ']';
    out.flags(flags);
    return out;
  }
  
  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &
  operator>>(std::basic_istream<char_t, traits_t> &in,
	     poisson_dist &g) {
    poisson_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed |
	     std::ios_base::left);
    in >> utility::ignore_spaces()
       >> utility::delim("[poisson ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }
  
}

#endif

