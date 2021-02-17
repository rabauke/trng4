// Copyright (c) 2000-2021, Heiko Bauke
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

#if !(defined TRNG_XOSHIRO256PLUS_HPP)

#define TRNG_XOSHIRO256PLUS_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/generate_canonical.hpp>
#include <trng/linear_algebra.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <ciso646>

namespace trng {

  class xoshiro256plus {
  public:
    // Uniform random number generator concept
    using result_type = uint64_t;
    TRNG_CUDA_ENABLE
    result_type operator()();

  private:
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = ~result_type(0);

  public:
    static constexpr result_type min() { return min_; }
    static constexpr result_type max() { return max_; }

    // Status class
    class status_type {
      result_type r[4]{result_type(1) << 63, 0, 0, 0};

    public:
      status_type() = default;
      explicit status_type(result_type r0, result_type r1, result_type r2, result_type r3)
          : r{r0, r1, r2, r3} {};

    private:
      struct eye {};
      explicit status_type(int i, eye) : status_type() {
        const result_type mask{result_type(1) << 63};
        if (i < 64)
          *this = status_type{mask >> i, 0, 0, 0};
        else if (i < 2 * 64)
          *this = status_type{0, mask >> (i - 64), 0, 0};
        else if (i < 3 * 64)
          *this = status_type{0, 0, mask >> (i - 2 * 64), 0};
        else
          *this = status_type{0, 0, 0, mask >> (i - 3 * 64)};
      }

    public:
      friend class xoshiro256plus;

      // Equality comparable concept
      friend bool operator==(const status_type &, const status_type &);
      friend bool operator!=(const status_type &, const status_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const status_type &S) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << S.r[0] << ' ' << S.r[1] << ' ' << S.r[2] << ' ' << S.r[3] << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, status_type &S) {
        status_type S_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> S_new.r[0] >> utility::delim(' ') >> S_new.r[1] >>
            utility::delim(' ') >> S_new.r[2] >> utility::delim(' ') >> S_new.r[3] >>
            utility::delim(')');
        if (in)
          S = S_new;
        in.flags(flags);
        return in;
      }

      vector<GF2, 256> to_vector() const {
        vector<GF2, 256> res;
        const result_type mask{result_type(1) << 63};
        for (int i{0}; i < 64; ++i)
          res(i) = GF2(((r[0] << i) & mask) > 0);
        for (int i{0}; i < 64; ++i)
          res(i + 64) = GF2(((r[1] << i) & mask) > 0);
        for (int i{0}; i < 64; ++i)
          res(i + 2 * 64) = GF2(((r[2] << i) & mask) > 0);
        for (int i{0}; i < 64; ++i)
          res(i + 3 * 64) = GF2(((r[3] << i) & mask) > 0);
        return res;
      }
    };

    // Random number engine concept
    explicit xoshiro256plus();
    explicit xoshiro256plus(unsigned long);
    explicit xoshiro256plus(result_type s0, result_type s1, result_type s2, result_type s3);

    template<typename gen>
    explicit xoshiro256plus(gen &g) {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      for (int k{0}; k < 4; ++k) {
        result_type r{0};
        for (int i{0}; i < 2; ++i) {
          r <<= 32u;
          r += g();
        }
        S.r[k] = r;
      }
      if (S.r[0] == 0 and S.r[1] == 0 and S.r[2] == 0 and S.r[3] == 0)
        S.r[0] |= result_type(1) << 63u;
    }
    void seed(result_type s0, result_type s1, result_type s2, result_type s3);

    // Equality comparable concept
    friend bool operator==(const xoshiro256plus &, const xoshiro256plus &);
    friend bool operator!=(const xoshiro256plus &, const xoshiro256plus &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const xoshiro256plus &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << xoshiro256plus::name() << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, xoshiro256plus &R) {
      xoshiro256plus::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(xoshiro256plus::name()) >>
          utility::delim(' ') >> S_new >> utility::delim(']');
      if (in) {
        R.S = S_new;
      }
      in.flags(flags);
      return in;
    }

    // Parallel random number generator concept
    //    TRNG_CUDA_ENABLE
    //    void split(unsigned int, unsigned int);
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
    status_type S;
    static const char *const name_str;

    TRNG_CUDA_ENABLE
    void step();
  };

  // Inline and template methods

  TRNG_CUDA_ENABLE
  inline void xoshiro256plus::step() {
    const result_type t{S.r[1] << 17};
    S.r[2] ^= S.r[0];
    S.r[3] ^= S.r[1];
    S.r[1] ^= S.r[2];
    S.r[0] ^= S.r[3];
    S.r[2] ^= t;
    S.r[3] = (S.r[3] << 45) | (S.r[3] >> (64 - 45));
  }

  TRNG_CUDA_ENABLE
  inline xoshiro256plus::result_type xoshiro256plus::operator()() {
    step();
    return S.r[0] + S.r[3];
  }

  TRNG_CUDA_ENABLE
  inline long xoshiro256plus::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, xoshiro256plus>(*this) * x);
  }

  // Parallel random number generator concept

  TRNG_CUDA_ENABLE
  inline void xoshiro256plus::jump2(unsigned int s) {
    matrix<GF2, 256> M;
    for (int i{0}; i < 256; ++i) {
      xoshiro256plus R;
      R.S = status_type(i, status_type::eye{});
      R.step();
      const vector<GF2, 256> v{R.S.to_vector()};
      for (int j{0}; j < 256; ++j)
        M(j, i) = v(j);
    }

    for (unsigned int i{0}; i < s; ++i)
      M = power(M, 2);
    const vector<GF2, 256> v0{S.to_vector()};
    const vector<GF2, 256> v1{M * v0};

    for (int j{0}; j < 4; ++j) {
      result_type r_j{0};
      for (int i{0}; i < 64; ++i) {
        r_j <<= 1;
        if (static_cast<bool>(v1(i + 64 * j)))
          r_j |= 1;
      }
      S.r[j] = r_j;
    }
  }

  TRNG_CUDA_ENABLE
  inline void xoshiro256plus::jump(unsigned long long s) {
    if (s < 16) {
      for (unsigned int i{0}; i < s; ++i)
        step();
    } else {
      matrix<GF2, 256> M;
      for (int i{0}; i < 256; ++i) {
        xoshiro256plus R;
        R.S = status_type(i, status_type::eye{});
        R.step();
        const vector<GF2, 256> v{R.S.to_vector()};
        for (int j{0}; j < 256; ++j)
          M(j, i) = v(j);
      }

      const vector<GF2, 256> v0{S.to_vector()};
      const vector<GF2, 256> v1{power(M, s) * v0};

      for (int j{0}; j < 4; ++j) {
        result_type r_j{0};
        for (int i{0}; i < 64; ++i) {
          r_j <<= 1;
          if (static_cast<bool>(v1(i + 64 * j)))
            r_j |= 1;
        }
        S.r[j] = r_j;
      }
    }
  }

  TRNG_CUDA_ENABLE
  inline void xoshiro256plus::discard(unsigned long long s) { jump(s); }

}  // namespace trng

#endif
