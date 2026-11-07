// Copyright (c) 2000-2026, Heiko Bauke
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

#if !(defined TRNG_MWC192X_HPP)

#define TRNG_MWC192X_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/generate_canonical.hpp>
#include <trng/mwc_detail.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#if defined _MSC_VER && __cplusplus <= 201703
#include <ciso646>
#endif

namespace trng {

  // Scrambled multiply-with-carry generator with 64-bit output and 192 bits of state, see
  // S. Vigna, "Scrambled Multiply-with-Carry Pseudorandom Number Generators".
  class mwc192x {
  public:
    // Uniform random number generator concept
    using result_type = uint64_t;
    TRNG_CUDA_ENABLE
    result_type operator()();

  private:
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = ~result_type(0);
    static constexpr result_type mult = 0xff4e24dd418f5433u;  // multiplier A

  public:
    TRNG_CUDA_ENABLE
    static constexpr result_type min() { return min_; }
    TRNG_CUDA_ENABLE
    static constexpr result_type max() { return max_; }

    // Status class
    class status_type {
      result_type x{0}, y{0}, c{1};

    public:
      status_type() = default;
      explicit status_type(result_type x, result_type y, result_type c) : x{x}, y{y}, c{c} {}

      friend class mwc192x;

      // Equality comparable concept
      friend bool operator==(const status_type &, const status_type &);
      friend bool operator!=(const status_type &, const status_type &);

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const status_type &S) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << S.x << ' ' << S.y << ' ' << S.c << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, status_type &S) {
        status_type S_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> S_new.x >> utility::delim(' ') >> S_new.y >>
            utility::delim(' ') >> S_new.c >> utility::delim(')');
        if (in)
          S = S_new;
        in.flags(flags);
        return in;
      }
    };

    // Random number engine concept
    explicit mwc192x();
    explicit mwc192x(unsigned long);
    explicit mwc192x(result_type x, result_type y, result_type c);

    template<typename gen>
    explicit mwc192x(gen &g) {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      S.x = draw(g);
      S.y = draw(g);
      S.c = clamp_carry(draw(g));
    }
    void seed(result_type x, result_type y, result_type c);

    // Equality comparable concept
    friend bool operator==(const mwc192x &, const mwc192x &);
    friend bool operator!=(const mwc192x &, const mwc192x &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const mwc192x &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << mwc192x::name() << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, mwc192x &R) {
      mwc192x::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(mwc192x::name()) >> utility::delim(' ') >>
          S_new >> utility::delim(']');
      if (in and mwc192x::canonical_state(S_new))
        R.S = S_new;
      else if (in)
        in.setstate(std::ios_base::failbit);
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

    // The carry must satisfy 0 < c < A-1; map invalid values into [1, A-2].
    static result_type clamp_carry(result_type c) {
      return c > 0u and c < mult - 1u ? c : c % (mult - 2u) + 1u;
    }


    template<typename gen>
    static result_type draw(gen &g) {
      result_type r{0};
      for (int i{0}; i < 2; ++i) {
        r <<= 32u;
        r += g();
      }
      return r;
    }


    TRNG_CUDA_ENABLE
    void step();

    // The state x + y*2^64 + c*2^128 is an integer modulo m = A*2^128 - 1.  A step
    // multiplies it by 2^-64 mod m, which is A*2^64, because A*2^128 = m + 1.
    using limbs_type = mwc_detail::limbs<6>;


    TRNG_CUDA_ENABLE
    static uint32_t lo32(const result_type x) { return static_cast<uint32_t>(x); }

    TRNG_CUDA_ENABLE
    static uint32_t hi32(const result_type x) { return static_cast<uint32_t>(x >> 32u); }

    TRNG_CUDA_ENABLE
    static result_type join(const uint32_t lo, const uint32_t hi) {
      return static_cast<result_type>(lo) | (static_cast<result_type>(hi) << 32u);
    }


    TRNG_CUDA_ENABLE
    static limbs_type modulus() {
      return limbs_type{{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu, lo32(mult - 1u),
                         hi32(mult - 1u)}};
    }

    TRNG_CUDA_ENABLE
    static bool canonical_state(const status_type &S) {
      const limbs_type s{{lo32(S.x), hi32(S.x), lo32(S.y), hi32(S.y), lo32(S.c), hi32(S.c)}};
      return !mwc_detail::geq(s, modulus());
    }

    TRNG_CUDA_ENABLE
    static limbs_type step_multiplier() {
      return limbs_type{{0u, 0u, lo32(mult), hi32(mult), 0u, 0u}};
    }

    TRNG_CUDA_ENABLE
    void jump_mul(const limbs_type &);
  };

  // Inline and template methods

  TRNG_CUDA_ENABLE
  inline void mwc192x::step() {
    result_type hi{0}, lo{0};
    mwc_detail::mul_add(mult, S.x, S.c, hi, lo);
    S.x = S.y;
    S.y = lo;
    S.c = hi;
  }

  TRNG_CUDA_ENABLE
  inline mwc192x::result_type mwc192x::operator()() {
    const result_type result{S.y ^ S.c};
    step();
    return result;
  }

  TRNG_CUDA_ENABLE
  inline long mwc192x::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, mwc192x>(*this) * x);
  }

  // Parallel random number generator concept

  TRNG_CUDA_ENABLE
  inline void mwc192x::jump_mul(const limbs_type &r) {
    limbs_type s{{lo32(S.x), hi32(S.x), lo32(S.y), hi32(S.y), lo32(S.c), hi32(S.c)}};
    mwc_detail::mul_mod(s, r, modulus());
    S.x = join(s.w[0], s.w[1]);
    S.y = join(s.w[2], s.w[3]);
    S.c = join(s.w[4], s.w[5]);
  }

  TRNG_CUDA_ENABLE
  inline void mwc192x::jump2(const unsigned int s) {
    limbs_type r;
    mwc_detail::pow2_mod(r, step_multiplier(), s, modulus());
    jump_mul(r);
  }

  TRNG_CUDA_ENABLE
  inline void mwc192x::jump(const unsigned long long s) {
    if (s < 64) {  // stepping is faster for short distances
      for (unsigned long long i{0}; i < s; ++i)
        step();
    } else {
      limbs_type r;
      mwc_detail::pow_mod(r, step_multiplier(), s, modulus());
      jump_mul(r);
    }
  }

  TRNG_CUDA_ENABLE
  inline void mwc192x::discard(const unsigned long long s) {
    jump(s);
  }

}  // namespace trng

#endif
