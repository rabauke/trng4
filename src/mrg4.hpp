// Copyright (c) 2000-2010, Heiko Bauke
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

#if !(defined TRNG_MRG4_HPP)

#define TRNG_MRG4_HPP

#include <ostream>
#include <istream>
#include <stdexcept>
#include <trng/utility.hpp>

namespace trng {
  
  class mrg4;
  
  class mrg4 {
  public:
  
    // Uniform random number generator concept
    typedef long result_type;
    result_type operator()() const;
  private:
    static const result_type modulus=2147483647l;
  public:
    static const result_type min=0l;
    static const result_type max=modulus-1;

    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
      result_type a1, a2, a3, a4;
    public:
      parameter_type() :
        a1(0), a2(0), a3(0), a4(0) { };
      parameter_type(result_type a1, result_type a2,
                     result_type a3, result_type a4) :
        a1(a1), a2(a2), a3(a3), a4(a4) { };

      friend class mrg4;

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
            << P.a1 << ' ' << P.a2 << ' ' << P.a3 << ' ' << P.a4
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
           >> P_new.a3 >> utility::delim(' ')
           >> P_new.a4 >> utility::delim(')');
        if (in)
          P=P_new;
        in.flags(flags);
        return in;
      }

    };

    class status_type {
      result_type r1, r2, r3, r4;
    public:
      status_type() : r1(0), r2(1), r3(1), r4(1) { };
      status_type(result_type r1, result_type r2, 
		  result_type r3, result_type r4) :
        r1(r1), r2(r2), r3(r3), r4(r4) { };

      friend class mrg4;

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
            << S.r1 << ' ' << S.r2 << ' ' << S.r3 << ' ' << S.r4
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
           >> S_new.r3 >> utility::delim(' ')
           >> S_new.r4 >> utility::delim(')');
        if (in)
          S=S_new;
        in.flags(flags);
        return in;
      }

    };

    static const parameter_type LEcuyer1;
    static const parameter_type LEcuyer2;

    // Random number engine concept
    explicit mrg4(parameter_type=LEcuyer1);
    explicit mrg4(unsigned long, parameter_type=LEcuyer1);
    template<typename gen>
    explicit mrg4(gen &g, parameter_type P=LEcuyer1) : P(P), S() {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      result_type r1=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      result_type r2=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      result_type r3=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      result_type r4=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      S.r1=r1;
      S.r2=r2;
      S.r3=r3;
      S.r4=r4;
    }
    void seed(result_type, result_type, result_type, result_type);

    // Equality comparable concept
    friend bool operator==(const mrg4 &, const mrg4 &);
    friend bool operator!=(const mrg4 &, const mrg4 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &
    operator<<(std::basic_ostream<char_t, traits_t> &out, const mrg4 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed |
                std::ios_base::left);
      out << '[' << mrg4::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &
    operator>>(std::basic_istream<char_t, traits_t> &in, mrg4 &R) {
      mrg4::parameter_type P_new;
      mrg4::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed |
               std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
         >> utility::delim(mrg4::name()) >> utility::delim(' ')
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
    void split(unsigned int, unsigned int);
    void jump2(unsigned int);
    void jump(unsigned long long);

    // Other useful methods
    static const char * name();
    long operator()(long) const;

  private:
    parameter_type P;
    mutable status_type S;
    static const char * const name_str;
    
    void backward();
    void step() const;
  };
    
  // Inline and template methods
  
  inline void mrg4::step() const {
    unsigned long long t(static_cast<unsigned long long>(P.a1)*
			 static_cast<unsigned long long>(S.r1)+
			 static_cast<unsigned long long>(P.a2)*
			 static_cast<unsigned long long>(S.r2)+
			 static_cast<unsigned long long>(P.a3)*
			 static_cast<unsigned long long>(S.r3)+
			 static_cast<unsigned long long>(P.a4)*
			 static_cast<unsigned long long>(S.r4));
    S.r4=S.r3;  S.r3=S.r2;  S.r2=S.r1;  S.r1=utility::modulo<modulus, 4>(t);
  }

  inline mrg4::result_type mrg4::operator()() const {
    step();
    return S.r1;
  }
  
  inline long mrg4::operator()(long x) const {
    return static_cast<long>(utility::uniformco<double, mrg4>(*this)*x);
  }

}
  
#endif
