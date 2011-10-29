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

#if !(defined TRNG_YARN2_HPP)

#define TRNG_YARN2_HPP

#include <ostream>
#include <istream>
#include <stdexcept>
#include <trng/utility.hpp>

namespace trng {
  
  class yarn2;
  
  class yarn2 {
  public:
  
    // Uniform random number generator concept
    typedef long result_type;
    result_type operator()() const;
  private:
    static const result_type modulus;
  public:
    static const result_type min;
    static const result_type max;
    
    // Parameter and status classes
    class parameter_type;
    class status_type;
    
    class parameter_type {
      result_type a1, a2;
      utility::power<2147483647l> g;
    public:
      parameter_type() :
        a1(0), a2(0), g() { };
      parameter_type(result_type a1, result_type a2,
		     result_type g) :
        a1(a1), a2(a2), g(g) { };

      friend class yarn2;

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
            << P.a1 << ' ' << P.a2 << ' ' << P.g
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
           >> P_new.g  >> utility::delim(')');
        if (in)
          P=P_new;
        in.flags(flags);
        return in;
      }

    };

   class status_type {
      result_type r1, r2;
    public:
      status_type() : r1(0), r2(1) { };
      status_type(result_type r1, result_type r2) :
        r1(r1), r2(r2) { };

      friend class yarn2;

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
            << S.r1 << ' ' << S.r2
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
           >> S_new.r2 >> utility::delim(')');
        if (in)
          S=S_new;
        in.flags(flags);
        return in;
      }

    };
  
    static const parameter_type LEcuyer1;
    static const parameter_type LEcuyer2;

    // Random number engine concept
    explicit yarn2(parameter_type=LEcuyer1);
    explicit yarn2(unsigned long, parameter_type=LEcuyer1);
    template<typename gen>
    explicit yarn2(gen &g, parameter_type P=LEcuyer1) : P(P), S() {
      result_type r1=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      result_type r2=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      S.r1=r1;
      S.r2=r2;
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      result_type r1=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      result_type r2=static_cast<unsigned long>(g())%
        static_cast<unsigned long>(modulus);
      S.r1=r1;
      S.r2=r2;
    }
    void seed(result_type, result_type);

    // Equality comparable concept
    friend bool operator==(const yarn2 &, const yarn2 &);
    friend bool operator!=(const yarn2 &, const yarn2 &);
    
    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &
    operator<<(std::basic_ostream<char_t, traits_t> &out, const yarn2 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed |
                std::ios_base::left);
      out << '[' << yarn2::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &
    operator>>(std::basic_istream<char_t, traits_t> &in, yarn2 &R) {
      yarn2::parameter_type P_new;
      yarn2::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed |
               std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[')
         >> utility::delim(yarn2::name()) >> utility::delim(' ')
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
    parameter_type P;
    mutable status_type S;
    static const char * const name_str;
    
    void backward();
    void step() const;
  };
    
  // Inline and template methods

  inline void yarn2::step() const {
    unsigned long long t(static_cast<unsigned long long>(P.a1)*
			 static_cast<unsigned long long>(S.r1)+
			 static_cast<unsigned long long>(P.a2)*
			 static_cast<unsigned long long>(S.r2));
    t=(t&0x7fffffffull)+(t>>31);
    if (t>=2ull*2147483647ull)
      t-=2ull*2147483647ull;
    if (t>=2147483647ull)
      t-=2147483647ull;
    S.r2=S.r1;  S.r1=t;
  }

  inline yarn2::result_type yarn2::operator()() const {
    step();
    return (S.r1==0) ? 0 : P.g(S.r1);
  }

  inline long yarn2::operator()(long x) const {
    return static_cast<long>(utility::uniformco(*this)*x);
  }
  
//   inline bool yarn2::boolean() const {
//     return operator()()<=modulus/2;
//   }

//   inline bool yarn2::boolean(double p) const {
//     return operator()()<modulus*p;
//   }
  
//   inline double yarn2::uniformco() const {
//     return static_cast<double>(operator()())/static_cast<double>(modulus);
//   }

//   inline double yarn2::uniformco(double a, double b) const {
//     return uniformco()*(b-a)+a;
//   }
  
//   inline double yarn2::uniformoc() const {
//     return (static_cast<double>(operator()())+1.0)/static_cast<double>(modulus);
//   }

//   inline double yarn2::uniformoc(double a, double b) const {
//     return uniformoc()*(b-a)+a;
//   }
  
//   inline double yarn2::uniformoo() const {
//     return (static_cast<double>(operator()())+1.0)/(static_cast<double>(modulus)+1.0);
//   }

//   inline double yarn2::uniformoo(double a, double b) const {
//     return uniformoo()*(b-a)+a;
//   }

//   inline double yarn2::uniformcc() const {
//     return static_cast<double>(operator()())/(static_cast<double>(modulus)-1.0);
//   }

//   inline double yarn2::uniformcc(double a, double b) const {
//     return uniformcc()*(b-a)+a;
//   }
    
}
  
#endif
