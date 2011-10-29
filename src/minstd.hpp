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

#if !(defined TRNG_MINSTD_HPP)

#define TRNG_MINSTD_HPP

#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <trng/utility.hpp>

namespace trng {
  
  class minstd;
  
  class minstd {
  public:

    // Uniform random number generator concept
    typedef unsigned long result_type;
    result_type operator()() const;  
    static const result_type min;
    static const result_type max;

    // Parameter and status classes
    class status_type;
    
    class status_type {
      result_type r;
    public:
      status_type() : r(1) { };
      explicit status_type(result_type r) : r(r) { };
      
      friend class minstd;

      // Equality comparable concept
      friend bool operator==(const status_type &, const status_type &);
      friend bool operator!=(const status_type &, const status_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, 
		 const status_type &S) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed | 
		  std::ios_base::left);
	out << '(' 
	    << S.r 
	    << ')';
	out.flags(flags);
	return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> &in, 
		 status_type &S) {
	status_type S_new;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed | 
		 std::ios_base::left);
	in >> utility::delim('(')
	   >> S_new.r >> utility::delim(')');
	if (in)
	  S=S_new;
	in.flags(flags);
	return in;
      }

    };
    
    // Random number engine concept
    explicit minstd(); 
    explicit minstd(unsigned long); 
    
    template<typename gen>
    explicit minstd(gen &g) : S() {
      seed(g);
    }
    
    void seed(); 
    void seed(result_type); 
    template<typename gen>
    void seed(gen &g) {
      do {
	S.r=g()%2147483647l;
	if (S.r<0)
	  S.r+=2147483647l;
      } while (S.r==0);
    }
    
    // Equality comparable concept
    friend bool operator==(const minstd &, const minstd &);
    friend bool operator!=(const minstd &, const minstd &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const minstd &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | 
		std::ios_base::left);
      out << '[' << minstd::name() << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, minstd &R) {
      minstd::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | 
	       std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
	 >> utility::delim(minstd::name()) >> utility::delim(' ')
	 >> S_new >> utility::delim(']');
      if (in) 
	R.S=S_new;
      in.flags(flags);
      return in;
    }
        
    // Other useful methods
    static const char * name();
    long operator()(long) const;
//     bool boolean() const;
//     bool boolean(double) const;
//     double uniformco() const;
//     double uniformco(double, double) const;
//     double uniformoc() const;
//     double uniformoc(double, double) const;
//     double uniformoo() const;
//     double uniformoo(double, double) const;
//     double uniformcc() const;
//     double uniformcc(double, double) const;
    
  private:
    mutable status_type S;
    static const char * const name_str;
    
    void step() const;
  };
  
  // Inline and template methods
  
  inline void minstd::step() const {
    unsigned long long t=S.r*16807;
     t=(t&0x7fffffffull)+(t>>31);
     if (t>=2147483647ull)
       t-=2147483647ull;
     S.r=static_cast<long>(t);
  }
  
  inline minstd::result_type minstd::operator()() const {
    step();
    return S.r;
  }

  inline long minstd::operator()(long x) const {
    return static_cast<long>(utility::uniformco(*this)*x);
  }

//   inline bool minstd::boolean() const {
//     step();
//     return S.r<9223372036854775808ull;
//   }
  
//   inline bool minstd::boolean(double p) const {
//     step();
//     return S.r<18446744073709551616.0*p;
//   }
  
//   inline double minstd::uniformco() const {
//     step();
//     double t(S.r/18446744073709551616.0);
//     return t<1.0 ? t : 1.0-math::numeric_limits<double>::epsilon();
//   }
  
//   inline double minstd::uniformco(double a, double b) const {
//     return uniformco()*(b-a)+a;
//     }
  
//   inline double minstd::uniformoc() const {
//     step();
//     double t(S.r/18446744073709551616.0);
//     return  t>0.0 ? t : math::numeric_limits<double>::epsilon();
//     }
  
//   inline double minstd::uniformoc(double a, double b) const {
//     return uniformoc()*(b-a)+a;
//   }
  
//   inline double minstd::uniformoo() const {
//     step();
//     double t(S.r/18446744073709551616.0+
// 	     math::numeric_limits<double>::epsilon());
//     return t<1.0 ? t : 1.0-math::numeric_limits<double>::epsilon();
//   }
  
//   inline double minstd::uniformoo(double a, double b) const {
//     return uniformoo()*(b-a)+a;
//   }
  
//   inline double minstd::uniformcc() const {
//     step();
//     return S.r/18446744073709551615.0;
//   }
  
//   inline double minstd::uniformcc(double a, double b) const {
//     return uniformcc()*(b-a)+a;
//   }

}

#endif
