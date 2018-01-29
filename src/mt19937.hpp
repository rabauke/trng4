// Copyright (c) 2000-2018, Heiko Bauke
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.  
// 
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.  
// 
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.


// This is a 32-bit version of Mersenne Twister pseudorandom number
// generator.
//
// References:
// M. Matsumoto and T. Nishimura,
//   ``Mersenne Twister: a 623-dimensionally equidistributed
//     uniform pseudorandom number generator''
//   ACM Transactions on Modeling and 
//   Computer Simulation 8. (Jan. 1998) 3--30.

#if !(defined TRNG_MT19937_HPP)

#define TRNG_MT19937_HPP

#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <trng/int_types.hpp>
#include <trng/utility.hpp>
#include <trng/generate_canonical.hpp>

namespace trng {
  
  class mt19937;
  
  class mt19937 {
  public:

    // Uniform random number generator concept
    typedef uint32_t result_type;
    result_type operator()();  
  private:
    static const result_type min_=0;
    static const result_type max_=4294967295u;
  public:
#if __cplusplus>=201103L
    static constexpr result_type min() {  return min_;  }
    static constexpr result_type max() {  return max_;  }
#else
    static const result_type min=min_;
    static const result_type max=max_;
#endif
  private:
    static const int N=624;
    static const int M=397;
    static const result_type UM=0x80000000u; // most significant bit
    static const result_type LM=0x7FFFFFFFu; // least significant 31 bits 
  public:

    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
    public:
      parameter_type() { };

      friend class mt19937;

      // Equality comparable concept
      friend bool operator==(const parameter_type &, const parameter_type &);
      friend bool operator!=(const parameter_type &, const parameter_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, 
		 const parameter_type &) {
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
      static const int N=624;
      int mti;
      result_type mt[N];
    public:
      status_type() : mti(0) { 
	for (int i=0; i<N; ++i)
	  mt[i]=0;
      }
      
      friend class mt19937;

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
    mt19937(); 
    explicit mt19937(unsigned long); 
    
    template<typename gen>
    explicit mt19937(gen &g) : P(), S() {
      seed(g);
    }
    
    void seed(); 
    template<typename gen>
    void seed(gen &g) {
      result_type r=g();
      seed(r);
    }
    void seed(unsigned long);

    void discard(unsigned long long);

    // Equality comparable concept
    friend bool operator==(const mt19937 &, const mt19937 &);
    friend bool operator!=(const mt19937 &, const mt19937 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const mt19937 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | 
		std::ios_base::left);
      out << '[' << mt19937::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, mt19937 &R) {
      mt19937::parameter_type P_new;
      mt19937::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | 
	       std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
	 >> utility::delim(mt19937::name()) >> utility::delim(' ')
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
    long operator()(long);

  private:
    parameter_type P;
    status_type S;
    static const char * const name_str;
    
  };
  
  // Inline and template methods
  
  inline mt19937::result_type mt19937::operator()() {
    result_type x;
    const result_type mag01[2]={ 0u, 0x9908b0dfu };
    if (S.mti>=N) { // generate N words at one time 
        int i;
        for (i=0; i<N-M; ++i) {
	  x=(S.mt[i]&mt19937::UM)|(S.mt[i+1]&mt19937::LM);
	  S.mt[i]=S.mt[i+M] ^ (x >> 1) ^ mag01[x & 0x1u];
        }
        for (; i<N-1; ++i) {
	  x=(S.mt[i]&mt19937::UM)|(S.mt[i+1]&LM);
	  S.mt[i]=S.mt[i+(M-N)] ^ (x >> 1) ^ mag01[x & 0x1u];
        }
        x=(S.mt[N-1]&mt19937::UM)|(S.mt[0]&mt19937::LM);
        S.mt[N-1]=S.mt[M-1] ^ (x >> 1) ^ mag01[x & 0x1u];
	S.mti=0;
    }
    x=S.mt[S.mti++];
    x^=(x >> 11);
    x^=(x << 7) & 0x9d2c5680u;
    x^=(x << 15) & 0xefc60000u;
    x^=(x >> 18);
    return x;
  }

  inline void mt19937::discard(unsigned long long n) {
    for (unsigned long long i(0); i<n; ++i)
      this->operator()();
  }

  inline long mt19937::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, mt19937>(*this)*x);
  }

}

#endif
