// Copyright (c) 2000-2020, Heiko Bauke
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
#include <trng/int_types.hpp>
#include <trng/utility.hpp>
#include <trng/generate_canonical.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <ciso646>

namespace trng {

  class mt19937_64 {
  public:
    // Uniform random number generator concept
    using result_type = uint64_t;
    result_type operator()();

  private:
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = 18446744073709551615u;

  public:
    static constexpr result_type min() { return min_; }
    static constexpr result_type max() { return max_; }

  private:
    static constexpr int N = 312;
    static constexpr int M = 156;
    static constexpr result_type UM = 0xFFFFFFFF80000000u;  // most significant 33 bits
    static constexpr result_type LM = 0x7FFFFFFFu;          // least significant 31 bits

  public:
    // Parameter and status classes
    class parameter_type {
    public:
      parameter_type() = default;

      friend class mt19937_64;

      // Equality comparable concept
      friend bool operator==(const parameter_type &, const parameter_type &);
      friend bool operator!=(const parameter_type &, const parameter_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const parameter_type &) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, parameter_type &P) {
        parameter_type P_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> utility::delim(')');
        if (in)
          P = P_new;
        in.flags(flags);
        return in;
      }
    };

    class status_type {
      static constexpr int N = 312;
      int mti{0};
      result_type mt[N]{};

    public:
      status_type() = default;

      friend class mt19937_64;

      // Equality comparable concept
      friend bool operator==(const status_type &, const status_type &);
      friend bool operator!=(const status_type &, const status_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const status_type &S) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << S.mti << ' ' << utility::make_io_range(S.mt, S.mt + N, " ") << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, status_type &S) {
        status_type S_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> S_new.mti >> utility::delim(' ') >>
            utility::make_io_range(S_new.mt, S_new.mt + N, " ") >> utility::delim(')');
        if (in)
          S = S_new;
        in.flags(flags);
        return in;
      }
    };

    // Random number engine concept
    mt19937_64();
    explicit mt19937_64(unsigned long);

    template<typename gen>
    explicit mt19937_64(gen &g) {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      result_type r{0};
      for (int i{0}; i < 2; ++i) {
        r <<= 32u;
        r += g();
      }
      seed(r);
    }

    void discard(unsigned long long);

    // Equality comparable concept
    friend bool operator==(const mt19937_64 &, const mt19937_64 &);
    friend bool operator!=(const mt19937_64 &, const mt19937_64 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const mt19937_64 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << mt19937_64::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, mt19937_64 &R) {
      mt19937_64::parameter_type P_new;
      mt19937_64::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(mt19937_64::name()) >> utility::delim(' ') >>
          P_new >> utility::delim(' ') >> S_new >> utility::delim(']');
      if (in) {
        R.P = P_new;
        R.S = S_new;
      }
      in.flags(flags);
      return in;
    }

    // Other useful methods
    static const char *name();
    long operator()(long);

  private:
    parameter_type P;
    status_type S;
    static const char *const name_str;
  };

  // Inline and template methods

  inline mt19937_64::result_type mt19937_64::operator()() {
    const result_type mag01[2]{0u, 0xB5026F5AA96619E9u};
    if (S.mti >= mt19937_64::N) {  // generate N words at one time
      int i{0};
      for (; i < mt19937_64::N - mt19937_64::M; ++i) {
        const result_type x{(S.mt[i] & mt19937_64::UM) | (S.mt[i + 1] & mt19937_64::LM)};
        S.mt[i] = S.mt[i + mt19937_64::M] ^ (x >> 1u) ^ mag01[(int)(x & 1u)];
      }
      for (; i < mt19937_64::N - 1; ++i) {
        const result_type x{(S.mt[i] & mt19937_64::UM) | (S.mt[i + 1] & mt19937_64::LM)};
        S.mt[i] = S.mt[i + (mt19937_64::M - mt19937_64::N)] ^ (x >> 1u) ^ mag01[(int)(x & 1u)];
      }
      const result_type x{(S.mt[mt19937_64::N - 1] & UM) | (S.mt[0] & LM)};
      S.mt[N - 1] = S.mt[mt19937_64::M - 1] ^ (x >> 1u) ^ mag01[(int)(x & 1u)];
      S.mti = 0;
    }
    result_type x{S.mt[S.mti++]};
    x ^= (x >> 29u) & 0x5555555555555555u;
    x ^= (x << 17u) & 0x71D67FFFEDA60000u;
    x ^= (x << 37u) & 0xFFF7EEE000000000u;
    x ^= (x >> 43u);
    return x;
  }

  inline void mt19937_64::discard(unsigned long long n) {
    for (unsigned long long i{0}; i < n; ++i)
      this->operator()();
  }

  inline long mt19937_64::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, mt19937_64>(*this) * x);
  }

}  // namespace trng

#endif
