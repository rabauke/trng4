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

#if !(defined TRNG_MRG3S_HPP)

#define TRNG_MRG3S_HPP

#include <ostream>
#include <istream>
#include <stdexcept>
#include <trng/cuda.hpp>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/int_math.hpp>

namespace trng {
  
  class mrg3s;
  
  class mrg3s {
  public:
    // Uniform random number generator concept
    typedef int32_t result_type;
    TRNG_CUDA_ENABLE
    result_type operator()();
  private:
    static const result_type modulus=2147462579;  // 2^31-21069
    static const result_type min_=0;
    static const result_type max_=modulus-1;
  public:
#if __cplusplus>=201103L
    static constexpr result_type min() {  return min_;  }
    static constexpr result_type max() {  return max_;  }
#else
    static const result_type min=min_;
    static const result_type max=max_;
#endif
    
    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
      result_type a1, a2, a3;
    public:
      parameter_type() :
	a1(0), a2(0), a3(0) { };
      parameter_type(result_type a1, result_type a2, 
		     result_type a3) :
	a1(a1), a2(a2), a3(a3) { };
      
      friend class mrg3s;
      
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
            << P.a1 << ' ' << P.a2 << ' ' << P.a3
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
           >> P_new.a1 >> utility::delim(' ')
           >> P_new.a2 >> utility::delim(' ')
           >> P_new.a3 >> utility::delim(')');
        if (in)
          P=P_new;
        in.flags(flags);
        return in;
      }
      
    };

    class status_type {
      result_type r1, r2, r3;
    public:
      status_type() : r1(0), r2(1), r3(1) { };
      status_type(result_type r1, result_type r2, 
		  result_type r3) : 
	r1(r1), r2(r2), r3(r3) { };
      
      friend class mrg3s;
      
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
            << S.r1 << ' ' << S.r2 << ' ' << S.r3
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
           >> S_new.r1 >> utility::delim(' ')
           >> S_new.r2 >> utility::delim(' ')
           >> S_new.r3 >> utility::delim(')');
        if (in)
          S=S_new;
        in.flags(flags);
        return in;
      }

    };
      
    static const parameter_type trng0;
    static const parameter_type trng1;
    
    // Random number engine concept
    explicit mrg3s(parameter_type=trng0);
    explicit mrg3s(unsigned long, parameter_type=trng0);
    template<typename gen>
    explicit mrg3s(gen &g, parameter_type P=trng0) : P(P), S() {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      result_type r1=static_cast<uint32_t>(g())%static_cast<uint32_t>(modulus);
      result_type r2=static_cast<uint32_t>(g())%static_cast<uint32_t>(modulus);
      result_type r3=static_cast<uint32_t>(g())%static_cast<uint32_t>(modulus);
      S.r1=r1;
      S.r2=r2;
      S.r3=r3;
    }
    void seed(result_type, result_type, result_type);

    // Equality comparable concept
    friend bool operator==(const mrg3s &, const mrg3s &);
    friend bool operator!=(const mrg3s &, const mrg3s &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &
    operator<<(std::basic_ostream<char_t, traits_t> &out, const mrg3s &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed |
                std::ios_base::left);
      out << '[' << mrg3s::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &
    operator>>(std::basic_istream<char_t, traits_t> &in, mrg3s &R) {
      mrg3s::parameter_type P_new;
      mrg3s::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed |
               std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
         >> utility::delim(mrg3s::name()) >> utility::delim(' ')
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
  inline void mrg3s::step() {
    uint64_t t(static_cast<uint64_t>(P.a1)*static_cast<uint64_t>(S.r1)+
	       static_cast<uint64_t>(P.a2)*static_cast<uint64_t>(S.r2)+
	       static_cast<uint64_t>(P.a3)*static_cast<uint64_t>(S.r3));
    S.r3=S.r2;  S.r2=S.r1;  S.r1=int_math::modulo<modulus, 3>(t);
  }

  TRNG_CUDA_ENABLE
  inline mrg3s::result_type mrg3s::operator()() {
    step();
    return S.r1;
  }
  
  TRNG_CUDA_ENABLE
  inline long mrg3s::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, mrg3s>(*this)*x);
  }
    
  // Parallel random number generator concept
  TRNG_CUDA_ENABLE
  inline void mrg3s::split(unsigned int s, unsigned int n) {
#if !(defined __CUDA_ARCH__)
    if (s<1 or n>=s)
      utility::throw_this(std::invalid_argument("invalid argument for trng::mrg3s::split"));
#endif
    if (s>1) {
      jump(n+1);  int32_t q0=S.r1;
      jump(s);    int32_t q1=S.r1;
      jump(s);    int32_t q2=S.r1;
      jump(s);    int32_t q3=S.r1;
      jump(s);    int32_t q4=S.r1;
      jump(s);    int32_t q5=S.r1;
      int32_t a[3], b[9];
      a[0]=q3;  b[0]=q2;  b[1]=q1;  b[2]=q0;
      a[1]=q4;  b[3]=q3;  b[4]=q2;  b[5]=q1;
      a[2]=q5;  b[6]=q4;  b[7]=q3;  b[8]=q2;
      int_math::gauss<3>(b, a, modulus);
      P.a1=a[0];  P.a2=a[1];  P.a3=a[2];
      S.r1=q2;    S.r2=q1;    S.r3=q0;
      for (int i=0; i<3; ++i)
	backward();
    }
  }
  
  TRNG_CUDA_ENABLE
  inline void mrg3s::jump2(unsigned int s) {
    int32_t b[9], c[9], d[3], r[3];
    int32_t t1(P.a1), t2(P.a2), t3(P.a3);
    b[0]=P.a1;  b[1]=P.a2;  b[2]=P.a3;
    b[3]=1;     b[4]=0;     b[5]=0;
    b[6]=0;     b[7]=1;     b[8]=0;
    for (unsigned int i(0); i<s; ++i)
      if ((i&1)==0)
	int_math::matrix_mult<3>(b, b, c, modulus);
      else
	int_math::matrix_mult<3>(c, c, b, modulus);
    r[0]=S.r1;  r[1]=S.r2;  r[2]=S.r3;
    if ((s&1)==0)
      int_math::matrix_vec_mult<3>(b, r, d, modulus);
    else
      int_math::matrix_vec_mult<3>(c, r, d, modulus);
    S.r1=d[0];  S.r2=d[1];  S.r3=d[2];
    P.a1=t1;    P.a2=t2;    P.a3=t3;
  }

  TRNG_CUDA_ENABLE
  inline void mrg3s::jump(unsigned long long s) {
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
  inline void mrg3s::discard(unsigned long long n) {
    return jump(n);
  }

  TRNG_CUDA_ENABLE
  inline void mrg3s::backward() {
    result_type t;
    if (P.a3!=0) {
      t=S.r1;
      t-=static_cast<result_type>((static_cast<int64_t>(P.a1)*
				   static_cast<int64_t>(S.r2))%modulus);
      if (t<0)
	t+=modulus;
      t-=static_cast<result_type>((static_cast<int64_t>(P.a2)*
				   static_cast<int64_t>(S.r3))%modulus);
      if (t<0)
	t+=modulus;
      t=static_cast<result_type>((static_cast<int64_t>(t)*
				  static_cast<int64_t>(int_math::modulo_invers(P.a3, modulus)))%modulus);
    } else if (P.a2!=0) {
      t=S.r2;
      t-=static_cast<result_type>((static_cast<int64_t>(P.a1)*
				   static_cast<int64_t>(S.r3))%modulus);
      if (t<0)
	t+=modulus;
      t=static_cast<result_type>((static_cast<int64_t>(t)*
				  static_cast<int64_t>(int_math::modulo_invers(P.a2, modulus)))%modulus);
    } else if (P.a1!=0) {
      t=S.r3;
      t=static_cast<result_type>((static_cast<int64_t>(t)*
				  static_cast<int64_t>(int_math::modulo_invers(P.a1, modulus)))%modulus);
    } else
      t=0;
    S.r1=S.r2;  S.r2=S.r3;  S.r3=t;
  }
  
}

#endif
