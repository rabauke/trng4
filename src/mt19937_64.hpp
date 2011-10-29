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


// This is a 64-bit version of Mersenne Twister pseudorandom number
// generator.
//
// References:
// T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
//   ACM Transactions on Modeling and 
//   Computer Simulation 10. (2000) 348--357.
// M. Matsumoto and T. Nishimura,
//   ``Mersenne Twister: a 623-dimensionally equidistributed
//     uniform pseudorandom number generator''
//   ACM Transactions on Modeling and 
//   Computer Simulation 8. (Jan. 1998) 3--30.

#if !(defined TRNG_MT19937_64_HPP)

#define TRNG_MT19937_64_HPP

#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <trng/utility.hpp>
#include <trng/generate_canonical.hpp>

namespace trng {
  
  class mt19937_64;
  
  class mt19937_64 {
  public:

    // Uniform random number generator concept
    typedef unsigned long long result_type;
    result_type operator()() const;  
    static const result_type min=0;
    static const result_type max=18446744073709551615ull;
  private:
    static const int N=312;
    static const int M=156;
    static const result_type UM=0xFFFFFFFF80000000ull; // most significant 33 bits
    static const result_type LM=0x7FFFFFFFull;         // least significant 31 bits 
  public:

    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
    public:
      parameter_type() { };

      friend class mt19937_64;

      // Equality comparable concept
      friend bool operator==(const parameter_type &, const parameter_type &);
      friend bool operator!=(const parameter_type &, const parameter_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, 
		 const parameter_type &P) {
	std::ios_base::fmtflags flags(out.flags());
	out.flags(std::ios_base::dec | std::ios_base::fixed | 
		  std::ios_base::left);
	out << '(' 
	    << ')';
	out.flags(flags);
	return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> &in, 
		 parameter_type &P) {
	parameter_type P_new;
	std::ios_base::fmtflags flags(in.flags());
	in.flags(std::ios_base::dec | std::ios_base::fixed | 
		 std::ios_base::left);
	in >> utility::delim('(') >> utility::delim(')');
	if (in)
	  P=P_new;
	in.flags(flags);
	return in;
      }
      
    };
    
    class status_type {
      static const int N=312;
      int mti;
      result_type mt[N];
    public:
      status_type() : mti(0) { 
	for (int i=0; i<N; ++i)
	  mt[i]=0;
      }
      
      friend class mt19937_64;

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
	    << S.mti << ' ' << utility::make_io_range(S.mt, S.mt+N, " ")
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
	   >> S_new.mti >> utility::delim(' ') >> utility::make_io_range(S_new.mt, S_new.mt+N, " ")
	   >> utility::delim(')');
	if (in)
	  S=S_new;
	in.flags(flags);
	return in;
      }

    };
    
    // Random number engine concept
    mt19937_64(); 
    explicit mt19937_64(unsigned long); 
    
    template<typename gen>
    explicit mt19937_64(gen &g) : P(), S() {
      seed(g);
    }
    
    void seed(); 
    void seed(unsigned long); 
    template<typename gen>
    void seed(gen &g) {
      result_type r=0;
      for (int i=0; i<2; ++i) {
	r<<=32;
	r+=g();
      }
      seed(r);
    }
    void seed(result_type);
    
    // Equality comparable concept
    friend bool operator==(const mt19937_64 &, const mt19937_64 &);
    friend bool operator!=(const mt19937_64 &, const mt19937_64 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const mt19937_64 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | 
		std::ios_base::left);
      out << '[' << mt19937_64::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, mt19937_64 &R) {
      mt19937_64::parameter_type P_new;
      mt19937_64::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | 
	       std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
	 >> utility::delim(mt19937_64::name()) >> utility::delim(' ')
	 >> P_new >> utility::delim(' ')
	 >> S_new >> utility::delim(']');
      if (in) { 
	R.P=P_new;
	R.S=S_new;
      }
      in.flags(flags);
      return in;
    }
    
    // Other useful methods
    static const char * name();
    long operator()(long) const;

  private:
    parameter_type P;
    mutable status_type S;
    static const char * const name_str;
    
  };
  
  // Inline and template methods
    
  inline mt19937_64::result_type mt19937_64::operator()() const {
    result_type x;
    const result_type mag01[2]={ 0ull, 0xB5026F5AA96619E9ull };
    
    if (S.mti>=mt19937_64::N) { // generate N words at one time
      int i;
      for (i=0; i<mt19937_64::N-mt19937_64::M; ++i) {
	x=(S.mt[i]&mt19937_64::UM)|(S.mt[i+1]&mt19937_64::LM);
	S.mt[i]=S.mt[i+mt19937_64::M] ^ (x>>1) ^ mag01[(int)(x&1ull)];
      }
      for (; i<mt19937_64::N-1; ++i) {
	x=(S.mt[i]&mt19937_64::UM)|(S.mt[i+1]&mt19937_64::LM);
	S.mt[i]=S.mt[i+(mt19937_64::M-mt19937_64::N)] ^ (x>>1) ^ mag01[(int)(x&1ull)];
      }
      x=(S.mt[mt19937_64::N-1]&UM)|(S.mt[0]&LM);
      S.mt[N-1]=S.mt[mt19937_64::M-1] ^ (x>>1) ^ mag01[(int)(x&1ull)];
      S.mti=0;
    }
    x=S.mt[S.mti++];
    x^=(x >> 29) & 0x5555555555555555ull;
    x^=(x << 17) & 0x71D67FFFEDA60000ull;
    x^=(x << 37) & 0xFFF7EEE000000000ull;
    x^=(x >> 43);
    return x;
  }

  inline long mt19937_64::operator()(long x) const {
    return static_cast<long>(utility::uniformco<double, mt19937_64>(*this)*x);
  }

}

#endif
