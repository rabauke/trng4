// Copyright (c) 2000-2024, Heiko Bauke
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
//     with the distribution
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

#if !(defined TRNG_COUNT128_LCG_SHIFT_HPP)

#define TRNG_COUNT128_LCG_SHIFT_HPP

#include <trng/trng_export.hpp>
#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/uint128.hpp>
#include <trng/generate_canonical.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <ciso646>

namespace trng {

  class count128_lcg_shift {
  public:
    // Uniform random number generator concept
    using result_type = uint64_t;
    TRNG_CUDA_ENABLE
    result_type operator()();

  private:
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = ~result_type(0);

  public:
    TRNG_CUDA_ENABLE
    static constexpr result_type min() { return min_; }
    TRNG_CUDA_ENABLE
    static constexpr result_type max() { return max_; }

  public:
    // Parameter and status classes
    class parameter_type {
      uint128 increment{0};
      result_type a{0}, b{0};

    public:
      parameter_type() = default;
      explicit parameter_type(uint128 increment, result_type a, result_type b)
          : increment{increment}, a{a}, b{b} {}

      friend class count128_lcg_shift;

      // Equality comparable concept
      friend bool operator==(const parameter_type &, const parameter_type &);
      friend bool operator!=(const parameter_type &, const parameter_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const parameter_type &P) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << P.increment << ' ' << P.a << ' ' << P.b << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, parameter_type &P) {
        parameter_type P_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> P_new.increment >> utility::delim(' ') >> P_new.a >>
            utility::delim(' ') >> P_new.b >> utility::delim(')');
        if (in)
          P = P_new;
        in.flags(flags);
        return in;
      }
    };

    class status_type {
      uint128 r{0};

    public:
      status_type() = default;
      explicit status_type(uint128 r) : r{r} {}

      friend class count128_lcg_shift;

      // Equality comparable concept
      friend bool operator==(const status_type &, const status_type &);
      friend bool operator!=(const status_type &, const status_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const status_type &S) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << S.r << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, status_type &S) {
        status_type S_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> S_new.r >> utility::delim(')');
        if (in)
          S = S_new;
        in.flags(flags);
        return in;
      }
    };

    static TRNG4_EXPORT const parameter_type Default;
    static TRNG4_EXPORT const parameter_type LEcuyer1;
    static TRNG4_EXPORT const parameter_type LEcuyer2;
    static TRNG4_EXPORT const parameter_type LEcuyer3;

    // Random number engine concept
    explicit count128_lcg_shift(parameter_type = Default);
    explicit count128_lcg_shift(unsigned long, parameter_type = Default);
    explicit count128_lcg_shift(unsigned long long, parameter_type = Default);

    template<typename gen>
    explicit count128_lcg_shift(gen &g, parameter_type P = Default) : P{P} {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      uint128 r{0};
      for (int i{0}; i < 4; ++i) {
        r <<= 32u;
        r += g();
      }
      S.r = r;
    }
    void seed(unsigned long long);

    // Equality comparable concept
    friend bool operator==(const count128_lcg_shift &, const count128_lcg_shift &);
    friend bool operator!=(const count128_lcg_shift &, const count128_lcg_shift &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const count128_lcg_shift &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << count128_lcg_shift::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, count128_lcg_shift &R) {
      count128_lcg_shift::parameter_type P_new;
      count128_lcg_shift::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(count128_lcg_shift::name()) >>
          utility::delim(' ') >> P_new >> utility::delim(' ') >> S_new >> utility::delim(']');
      if (in) {
        R.P = P_new;
        R.S = S_new;
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
    static const char *name();
    TRNG_CUDA_ENABLE
    long operator()(long);

  private:
    parameter_type P;
    status_type S;
    static const char *const name_str;

    TRNG_CUDA_ENABLE
    void backward();
    TRNG_CUDA_ENABLE
    void step();
  };

  // Inline and template methods

  TRNG_CUDA_ENABLE
  inline void count128_lcg_shift::step() {
    S.r += P.increment;
  }

  TRNG_CUDA_ENABLE
  inline count128_lcg_shift::result_type count128_lcg_shift::operator()() {
    step();
    const result_type t_hi{S.r.hi()};
    const result_type t_lo{S.r.lo()};
    result_type t{(t_lo ^ t_hi) * P.a + P.b};
    t ^= (t >> 23u);
    t ^= (t << 41u);
    t ^= (t >> 18u);
    return t;
  }

  TRNG_CUDA_ENABLE
  inline long count128_lcg_shift::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, count128_lcg_shift>(*this) * x);
  }

  // Parallel random number generator concept

  TRNG_CUDA_ENABLE
  inline void count128_lcg_shift::split(unsigned int s, unsigned int n) {
#if !(defined TRNG_CUDA)
    if (s < 1 or n >= s)
      utility::throw_this(
          std::invalid_argument("invalid argument for trng::count128_shift_lcg::split"));
#endif
    if (s > 1) {
      S.r += uint128{0, n} * P.increment;
      S.r += P.increment;
      P.increment *= uint128{0, s};
      S.r -= P.increment;
    }
  }

  TRNG_CUDA_ENABLE
  inline void count128_lcg_shift::jump2(unsigned int s) {
    S.r += (uint128{0, 1} << (s % 128)) * P.increment;
  }

  TRNG_CUDA_ENABLE
  inline void count128_lcg_shift::jump(unsigned long long s) {
    S.r += uint128{0, s} * P.increment;
  }

  TRNG_CUDA_ENABLE
  inline void count128_lcg_shift::discard(unsigned long long n) {
    jump(n);
  }

}  // namespace trng

#endif
