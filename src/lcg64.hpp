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

#if !(defined TRNG_LCG64_HPP)

#define TRNG_LCG64_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/generate_canonical.hpp>

namespace trng {
  
  class lcg64;
  
  class lcg64 {
  public:

    // Uniform random number generator concept
    typedef uint64_t result_type;
    TRNG_CUDA_ENABLE
    result_type operator()(); 
  private:
    static const result_type min_=0;
    static const result_type max_=18446744073709551615u;
  public:
#if __cplusplus>=201103L
    static constexpr result_type min() {  return min_;  }
    static constexpr result_type max() {  return max_;  }
#else
    static const result_type min=min_;
    static const result_type max=max_;
#endif
  private:
    // compute floor(log_2(x))
    TRNG_CUDA_ENABLE
    static unsigned int log2_floor(result_type x);
    // compute pow(x, n)
    TRNG_CUDA_ENABLE
    static result_type pow(result_type x, result_type n);
    // compute prod(1+a^(2^i), i=0..l-1)
    TRNG_CUDA_ENABLE
    static result_type g(unsigned int l, result_type a);
    // compute sum(a^i, i=0..s-1)
    TRNG_CUDA_ENABLE
    static result_type f(result_type s, result_type a);
  public:
    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
      result_type a, b;
    public:
      parameter_type() :
	a(0), b(0) { };
      parameter_type(result_type a, result_type b) :
	a(a), b(b) { };

      friend class lcg64;

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
	    << P.a << ' ' << P.b 
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
	in >> utility::delim('(')
	   >> P_new.a >> utility::delim(' ')
	   >> P_new.b >> utility::delim(')');
	if (in)
	  P=P_new;
	in.flags(flags);
	return in;
      }
      
    };
    
    class status_type {
      result_type r;
    public:
      status_type() : r(0) { };
      explicit status_type(result_type r) : r(r) { };
      
      friend class lcg64;

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
    
    static const parameter_type Default;
    static const parameter_type LEcuyer1;
    static const parameter_type LEcuyer2;
    static const parameter_type LEcuyer3;

    // Random number engine concept
    explicit lcg64(parameter_type=Default); 
    explicit lcg64(unsigned long, parameter_type=Default); 
    explicit lcg64(unsigned long long, parameter_type=Default); 
    
    template<typename gen>
    explicit lcg64(gen &g, parameter_type P=Default) : P(P), S() {
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
      S.r=r;
    }
    void seed(unsigned long long); 
    
    // Equality comparable concept
    friend bool operator==(const lcg64 &, const lcg64 &);
    friend bool operator!=(const lcg64 &, const lcg64 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const lcg64 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | 
		std::ios_base::left);
      out << '[' << lcg64::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, lcg64 &R) {
      lcg64::parameter_type P_new;
      lcg64::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | 
	       std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
	 >> utility::delim(lcg64::name()) >> utility::delim(' ')
	 >> P_new >> utility::delim(' ')
	 >> S_new >> utility::delim(']');
      if (in) { 
	R.P=P_new;
	R.S=S_new;
      }
      in.flags(flags);
      return in;
    }
    
    // Parallel random number generator concept
    TRNG_CUDA_ENABLE
    void split(unsigned int, unsigned int);
    TRNG_CUDA_ENABLE
    void jump2(unsigned int);
    TRNG_CUDA_ENABLE
    void jump(unsigned long long);
    TRNG_CUDA_ENABLE
    void discard(unsigned long long);
    
    // Other useful methods
    static const char * name();
    TRNG_CUDA_ENABLE
    long operator()(long);

  private:
    parameter_type P;
    status_type S;
    static const char * const name_str;
    
    TRNG_CUDA_ENABLE
    void backward();
    TRNG_CUDA_ENABLE
    void step();
  };
  
  // Inline and template methods
  
  TRNG_CUDA_ENABLE
  inline void lcg64::step() {
    S.r=(P.a*S.r+P.b);
  }
    
  TRNG_CUDA_ENABLE
  inline lcg64::result_type lcg64::operator()() {
    step();
    return S.r;
  }

  TRNG_CUDA_ENABLE
  inline long lcg64::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, lcg64>(*this)*x);
  }

  // compute floor(log_2(x))
  TRNG_CUDA_ENABLE
  inline unsigned int lcg64::log2_floor(lcg64::result_type x) {
    unsigned int y(0);
    while (x>0) {
      x>>=1;
      ++y;
    };
    --y;
    return y;
  }

  // compute pow(x, n)
  TRNG_CUDA_ENABLE
  inline lcg64::result_type lcg64::pow(lcg64::result_type x, lcg64::result_type n) {
    lcg64::result_type result=1;
    while (n>0) {
      if ((n&1)>0)
	result=result*x;
      x=x*x;
      n>>=1;
    }
    return result;
  }

  // compute prod(1+a^(2^i), i=0..l-1)
  TRNG_CUDA_ENABLE
  inline lcg64::result_type lcg64::g(unsigned int l, lcg64::result_type a) {
    lcg64::result_type p=a, res=1;
    for (unsigned i=0; i<l; ++i) {
      res*=1+p;
      p*=p;
    }
    return res;
  }
  
  // compute sum(a^i, i=0..s-1)
  TRNG_CUDA_ENABLE
  inline lcg64::result_type lcg64::f(lcg64::result_type s, lcg64::result_type a) {
    if (s==0)
      return 0;
    unsigned int e=log2_floor(s);
    lcg64::result_type y=0, p=a;
    for (unsigned int l=0; l<=e; ++l) {
      if (((1ull<<l) & s)>0) {
        y=g(l, a)+p*y;
      }
      p*=p;
    }
    return y;
  }

  // Parallel random number generator concept

  TRNG_CUDA_ENABLE
  inline void lcg64::split(unsigned int s, unsigned int n) {
#if !(defined __CUDA_ARCH__)
    if (s<1 or n>=s)
      utility::throw_this(std::invalid_argument("invalid argument for trng::lcg64::split"));
#endif
    if (s>1) {
      jump(n+1);
      P.b*=f(s, P.a);
      P.a=pow(P.a, s);
      backward();
    }
  }

  TRNG_CUDA_ENABLE
  inline void lcg64::jump2(unsigned int s) {
    S.r=S.r*pow(P.a, 1ull<<s)+f(1ull<<s, P.a)*P.b;
  }

  TRNG_CUDA_ENABLE
  inline void lcg64::jump(unsigned long long s) {
    if (s<16) {
      for (unsigned int i(0); i<s; ++i) 
	step();
    } else {
      unsigned int i(0);
      while (s>0) {
	if (s%2==1)
	  jump2(i);
	++i;
	s>>=1;
      }
    }
  }
  
  TRNG_CUDA_ENABLE
  inline void lcg64::discard(unsigned long long n) {
    return jump(n);
  }

  TRNG_CUDA_ENABLE
  inline void lcg64::backward() {
    for (unsigned int i(0); i<64; ++i)
      jump2(i);
  }

}

#endif
