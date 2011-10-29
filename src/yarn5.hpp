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

#if !(defined TRNG_YARN5_HPP)

#define TRNG_YARN5_HPP

#include <ostream>
#include <istream>
#include <stdexcept>
#include <trng/utility.hpp>

namespace trng {
  
  class yarn5;
  
  class yarn5 {
  public:
  
    // Uniform random number generator concept
    typedef long result_type;
    result_type operator()() const;
  private:
    static const result_type modulus=2147483647l;
  public:
    static const result_type min=0l;
    static const result_type max=2147483646l;

    // Parameter and status classes
    class parameter_type;
    class status_type;

    class parameter_type {
      result_type a1, a2, a3, a4, a5;
      utility::power<2147483647l> g;
    public:
      parameter_type() :
        a1(0), a2(0), a3(0), a4(0), a5(0), g() { };
      parameter_type(result_type a1, result_type a2,
                     result_type a3, result_type a4,
		     result_type a5,
                     result_type g) :
        a1(a1), a2(a2), a3(a3), a4(a4), a5(a5), g(g) { };

      friend class yarn5;

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
            << P.a1 << ' ' << P.a2 << ' ' << P.a3 << ' ' << P.a4 << ' ' << P.a5 << ' ' << P.g
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
           >> P_new.a4 >> utility::delim(' ')
           >> P_new.a5 >> utility::delim(' ')
           >> P_new.g  >> utility::delim(')');
        if (in)
          P=P_new;
        in.flags(flags);
        return in;
      }

    };

    class status_type {
      result_type r1, r2, r3, r4, r5;
    public:
      status_type() : r1(0), r2(1), r3(1), r4(1), r5(1) { };
      status_type(result_type r1, result_type r2, 
		  result_type r3, result_type r4, 
		  result_type r5) :
        r1(r1), r2(r2), r3(r3), r4(r4), r5(r5) { };
      
      friend class yarn5;

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
            << S.r1 << ' ' << S.r2 << ' ' << S.r3 << ' ' << S.r4 << ' ' << S.r5
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
           >> S_new.r4 >> utility::delim(' ')
           >> S_new.r5 >> utility::delim(')');
        if (in)
          S=S_new;
        in.flags(flags);
        return in;
      }

    };
     
    static const parameter_type LEcuyer1;
    
    // Random number engine concept
    explicit yarn5(parameter_type=LEcuyer1);
    explicit yarn5(unsigned long, parameter_type=LEcuyer1);
    template<typename gen>
    explicit yarn5(gen &g, parameter_type P=LEcuyer1) : P(P), S() {
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
      result_type r5=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      S.r1=r1;
      S.r2=r2;
      S.r3=r3;
      S.r4=r4;
      S.r5=r5;
    }
    void seed(result_type, result_type, result_type, result_type, result_type);

    // Equality comparable concept
    friend bool operator==(const yarn5 &, const yarn5 &);
    friend bool operator!=(const yarn5 &, const yarn5 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &
    operator<<(std::basic_ostream<char_t, traits_t> &out, const yarn5 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed |
                std::ios_base::left);
      out << '[' << yarn5::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &
    operator>>(std::basic_istream<char_t, traits_t> &in, yarn5 &R) {
      yarn5::parameter_type P_new;
      yarn5::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed |
               std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
         >> utility::delim(yarn5::name()) >> utility::delim(' ')
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

  inline void yarn5::step() const {
    unsigned long long t(static_cast<unsigned long long>(P.a1)*
			 static_cast<unsigned long long>(S.r1)+
			 static_cast<unsigned long long>(P.a2)*
			 static_cast<unsigned long long>(S.r2)+
			 static_cast<unsigned long long>(P.a3)*
			 static_cast<unsigned long long>(S.r3)+
			 static_cast<unsigned long long>(P.a4)*
			 static_cast<unsigned long long>(S.r4));
    if (t>=2ull*2147483647ull*2147483647ull)
      t-=2ull*2147483647ull*2147483647ull;
    t+=static_cast<unsigned long long>(P.a5)*
      static_cast<unsigned long long>(S.r5);
    t=(t&0x7fffffffull)+(t>>31);
    if (t>=4ull*2147483647ull)
      t-=4ull*2147483647ull;
    if (t>=2ull*2147483647ull)
      t-=2ull*2147483647ull;
    if (t>=2147483647ull)
      t-=2147483647ull;
    S.r5=S.r4;  S.r4=S.r3;  S.r3=S.r2;  S.r2=S.r1;  S.r1=t;
  }

  inline yarn5::result_type yarn5::operator()() const {
    step();
    return (S.r1==0) ? 0 : P.g(S.r1);
  }
  
  inline long yarn5::operator()(long x) const {
    return static_cast<long>(utility::uniformco(*this)*x);
  }

}
  
#endif
