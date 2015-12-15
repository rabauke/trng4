// Copyright (c) 2000-2015, Heiko Bauke
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

#if !(defined TRNG_MINSTD_HPP)

#define TRNG_MINSTD_HPP

#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>

namespace trng {
  
  class minstd;
  
  class minstd {
  public:

    // Uniform random number generator concept
    typedef uint32_t result_type;
    result_type operator()();
  private:
    static const result_type min_=1u;
    static const result_type max_=2147483646u;
  public:
#if __cplusplus>=201103L
    static constexpr result_type min() {  return min_;  }
    static constexpr result_type max() {  return max_;  }
#else
    static const result_type min=min_;
    static const result_type max=max_;
#endif
    
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
    minstd(); 
    explicit minstd(unsigned long); 
    
    template<typename gen>
    explicit minstd(gen &g) : S() {
      seed(g);
    }
    
    void seed(); 
    void seed(unsigned long); 
    template<typename gen>
    void seed(gen &g) {
      do {
	S.r=g()%2147483647;
      } while (S.r==0);
    }
    
    void discard(unsigned long long);
    
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
    long operator()(long);
    
  private:
    status_type S;
    static const char * const name_str;
    
    void step();
  };
  
  // Inline and template methods
  
  inline void minstd::step() {
    uint64_t t=S.r*16807;
     t=(t&0x7fffffffu)+(t>>31);
     if (t>=2147483647u)
       t-=2147483647u;
     S.r=static_cast<long>(t);
  }
  
  inline minstd::result_type minstd::operator()() {
    step();
    return S.r;
  }

  inline void minstd::discard(unsigned long long n) {
    n%=minstd::max_-minstd::min_;
    for (unsigned long long i(0); i<n; ++i)
      step();
  }
  
  inline long minstd::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, minstd>(*this)*x);
  }
  
}

#endif
