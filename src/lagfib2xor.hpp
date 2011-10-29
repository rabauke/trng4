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

#if !(defined TRNG_LAGFIB2XOR_HPP)

#define TRNG_LAGFIB2XOR_HPP

#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <sstream>
#include <trng/utility.hpp>
#include <trng/minstd.hpp>

namespace trng {

  template<typename integer_type,
	   unsigned int A, unsigned int B>
  class lagfib2xor;
  
  template<typename integer_type,
	   unsigned int A, unsigned int B>
  class lagfib2xor {
  public:
    
    // Uniform random number generator concept
    typedef integer_type result_type;
    result_type operator()() const {
      step();  
      return S.r[S.index];
    }
    static const result_type min=0;
    static const result_type max=~result_type(0);
    
    // Parameter and status classes
    class status_type;

    class status_type {
      result_type r[utility::ceil2<B>::result];
      unsigned int index;
    public:
      status_type() { 
	for (unsigned int i=0; i<utility::ceil2<B>::result; ++i)
	  r[i]=0;
	index=0;
      };
      
      friend class lagfib2xor;
      
      // Equality comparable concept
      friend bool operator==(const status_type &a, const status_type &b) {
	if (a.index!=b.index) 
	  return false;
	for (unsigned int i=0; i<utility::ceil2<B>::result; ++i)
	  if (a.r[i]!=b.r[i])
	    return false;
	return true;
      }
      friend bool operator!=(const status_type &a, const status_type &b) {
	return not (a==b);
      }

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, 
		 const status_type &S) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed | 
		  std::ios_base::left);
	out << '(' 
	    << S.index;
 	for (unsigned int i=0; i<utility::ceil2<B>::result; ++i)
	  out << ' ' << S.r[i];
	out << ')';
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
	   >> S_new.index;
 	for (unsigned int i=0; i<utility::ceil2<B>::result; ++i)
 	  in >> utility::delim(' ') >> S_new.r[i];
	in >> utility::delim(')');
	if (in)
	  S=S_new;
	in.flags(flags);
	return in;
      }
      
    };
    
    // Random number engine concept
    lagfib2xor() : S() {
      seed();
    }
    
    explicit lagfib2xor(unsigned long s) : S() {
      seed(s);
    }
    
    template<typename gen>
    explicit lagfib2xor(gen &g) : S() {
      seed(g);
    }
    
    void seed() {
      seed(0);
    }
    
    void seed(unsigned long s) {
      minstd R(s);
      seed(R);
    }
    
    template<typename gen>
    void seed(gen &g) {
      for (unsigned int i=0; i<B; ++i) {
        result_type r=0;
        for (unsigned int j=0; j<std::numeric_limits<result_type>::digits; ++j) {
          r<<=1;
	  if (g()-gen::min>gen::max/2)
            ++r;
        }
        S.r[i]=r;
      }
      S.index=B-1;
    }
    
    // Equality comparable concept
    friend bool operator==(const lagfib2xor &R1, const lagfib2xor &R2) {
      return R1.S==R2.S;
    }
      
    friend bool operator!=(const lagfib2xor &R1, const lagfib2xor &R2) {
      return not (R1==R2);
    }
    
    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const lagfib2xor &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | 
		std::ios_base::left);
      out << '[' << lagfib2xor::name() << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }
    
    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, lagfib2xor &R) {
      typename lagfib2xor::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | 
	       std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
	 >> utility::delim(lagfib2xor::name()) >> utility::delim(' ')
	 >> S_new >> utility::delim(']');
      if (in)
	R.S=S_new;
      in.flags(flags);
      return in;
    }
    
    // Other useful methods
    static const char * name() {
      static char name_c_str[64]={'\0'};
      if (name_c_str[0]=='\0') {
	std::stringstream name_str;
	name_str << "lagfib2xor_" << std::numeric_limits<result_type>::digits << '_' 
		 << A << '_' << B;
	int i=0;
	const char *p=name_str.str().c_str();
	while (p[i]!='\0' and i<63) {
	  name_c_str[i]=p[i];
	  ++i;
	}
	name_c_str[i]='\0';
      }
      return name_c_str;
    }
    long operator()(long x) const {
      return static_cast<long>(utility::uniformco<double, lagfib2xor>(*this)*x);
    }
    
  private:
    mutable status_type S;
    
    void step() const {
      S.index++;
      S.index&=utility::mask<B>::result;
      S.r[S.index]=
	S.r[(S.index-A)&utility::mask<B>::result] ^ 
	S.r[(S.index-B)&utility::mask<B>::result];
    }
  };
  
  typedef lagfib2xor<unsigned long,       103,   250> r250_ul;
  typedef lagfib2xor<unsigned long long,  103,   250> r250_ull;
  typedef lagfib2xor<unsigned long,       168,   521> lagfib2xor_521_ul;
  typedef lagfib2xor<unsigned long long,  168,   521> lagfib2xor_521_ull;
  typedef lagfib2xor<unsigned long,       273,   607> lagfib2xor_607_ul;
  typedef lagfib2xor<unsigned long long,  273,   607> lagfib2xor_607_ull;
  typedef lagfib2xor<unsigned long,       418,  1279> lagfib2xor_1279_ul;
  typedef lagfib2xor<unsigned long long,  418,  1279> lagfib2xor_1279_ull;
  typedef lagfib2xor<unsigned long,      1029,  2281> lagfib2xor_2281_ul;
  typedef lagfib2xor<unsigned long long, 1029,  2281> lagfib2xor_2281_ull;
  typedef lagfib2xor<unsigned long,       576,  3217> lagfib2xor_3217_ul;
  typedef lagfib2xor<unsigned long long,  576,  3217> lagfib2xor_3217_ull;
  typedef lagfib2xor<unsigned long,      2098,  4423> lagfib2xor_4423_ul;
  typedef lagfib2xor<unsigned long long, 2098,  4423> lagfib2xor_4423_ull;
  typedef lagfib2xor<unsigned long,      4187,  9689> lagfib2xor_9689_ul;
  typedef lagfib2xor<unsigned long long, 4187,  9689> lagfib2xor_9689_ull;
  typedef lagfib2xor<unsigned long,      9842, 19937> lagfib2xor_19937_ul;
  typedef lagfib2xor<unsigned long long, 9842, 19937> lagfib2xor_19937_ull;
		       
}

#endif
