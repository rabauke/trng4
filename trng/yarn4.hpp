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

#if !(defined TRNG_YARN4_HPP)

#define TRNG_YARN4_HPP

#include <trng/cuda.hpp>
#include <trng/utility.hpp>
#include <trng/int_types.hpp>
#include <trng/int_math.hpp>
#include <trng/mrg_parameter.hpp>
#include <trng/mrg_status.hpp>
#include <ostream>
#include <istream>
#include <stdexcept>
#include <ciso646>

namespace trng {

  class yarn4 {
  public:
    // Uniform random number generator concept
    using result_type = int32_t;
    TRNG_CUDA_ENABLE
    result_type operator()();

  private:
    static constexpr result_type modulus = 2147483647;
    static constexpr result_type gen = 123567893;
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = modulus - 1;
    static const int_math::power<yarn4::modulus, yarn4::gen> g;

  public:
    static constexpr result_type min() { return min_; }
    static constexpr result_type max() { return max_; }

    // Parameter and status classes
    using parameter_type = mrg_parameter<result_type, 4, yarn4>;
    using status_type = mrg_status<result_type, 4, yarn4>;

    static const parameter_type LEcuyer1;
    static const parameter_type LEcuyer2;

    // Random number engine concept
    explicit yarn4(parameter_type = LEcuyer1);
    explicit yarn4(unsigned long, parameter_type = LEcuyer1);
    template<typename gen>
    explicit yarn4(gen &g, parameter_type P = LEcuyer1) : P{P} {
      seed(g);
    }

    void seed();
    void seed(unsigned long);
    template<typename gen>
    void seed(gen &g) {
      const result_type r1{static_cast<uint32_t>(g()) % static_cast<uint32_t>(modulus)};
      const result_type r2{static_cast<uint32_t>(g()) % static_cast<uint32_t>(modulus)};
      const result_type r3{static_cast<uint32_t>(g()) % static_cast<uint32_t>(modulus)};
      const result_type r4{static_cast<uint32_t>(g()) % static_cast<uint32_t>(modulus)};
      S.r[0] = r1;
      S.r[1] = r2;
      S.r[2] = r3;
      S.r[3] = r4;
    }
    void seed(result_type, result_type, result_type, result_type);

    // Equality comparable concept
    friend bool operator==(const yarn4 &, const yarn4 &);
    friend bool operator!=(const yarn4 &, const yarn4 &);

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const yarn4 &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << yarn4::name() << ' ' << R.P << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, yarn4 &R) {
      yarn4::parameter_type P_new;
      yarn4::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(yarn4::name()) >> utility::delim(' ') >>
          P_new >> utility::delim(' ') >> S_new >> utility::delim(']');
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
  inline void yarn4::step() {
    const uint64_t t{static_cast<uint64_t>(P.a[0]) * static_cast<uint64_t>(S.r[0]) +
                     static_cast<uint64_t>(P.a[1]) * static_cast<uint64_t>(S.r[1]) +
                     static_cast<uint64_t>(P.a[2]) * static_cast<uint64_t>(S.r[2]) +
                     static_cast<uint64_t>(P.a[3]) * static_cast<uint64_t>(S.r[3])};
    S.r[3] = S.r[2];
    S.r[2] = S.r[1];
    S.r[1] = S.r[0];
    S.r[0] = int_math::modulo<modulus, 4>(t);
  }

  TRNG_CUDA_ENABLE
  inline yarn4::result_type yarn4::operator()() {
    step();
#if defined __CUDA_ARCH__
    if (S.r[0] == 0)
      return 0;
    yarn4::result_type n = S.r[0];
    int64_t p(1), t(gen);
    while (n > 0) {
      if ((n & 0x1) == 0x1)
        p = int_math::modulo<modulus, 1>(p * t);
      t = int_math::modulo<modulus, 1>(t * t);
      n /= 2;
    }
    return static_cast<yarn4::result_type>(p);
#else
    return S.r[0] == 0 ? 0 : g(S.r[0]);
#endif
  }

  TRNG_CUDA_ENABLE
  inline long yarn4::operator()(long x) {
    return static_cast<long>(utility::uniformco<double, yarn4>(*this) * x);
  }

  // Parallel random number generator concept
  TRNG_CUDA_ENABLE
  inline void yarn4::split(unsigned int s, unsigned int n) {
#if !(defined __CUDA_ARCH__)
    if (s < 1 or n >= s)
      utility::throw_this(std::invalid_argument("invalid argument for trng::yarn4::split"));
#endif
    if (s > 1) {
      jump(n + 1);
      const int32_t q0{S.r[0]};
      jump(s);
      const int32_t q1{S.r[0]};
      jump(s);
      const int32_t q2{S.r[0]};
      jump(s);
      const int32_t q3{S.r[0]};
      jump(s);
      const int32_t q4{S.r[0]};
      jump(s);
      const int32_t q5{S.r[0]};
      jump(s);
      const int32_t q6{S.r[0]};
      jump(s);
      const int32_t q7{S.r[0]};
      int32_t a[4], b[16];
      a[0] = q4;
      b[0] = q3;
      b[1] = q2;
      b[2] = q1;
      b[3] = q0;
      a[1] = q5;
      b[4] = q4;
      b[5] = q3;
      b[6] = q2;
      b[7] = q1;
      a[2] = q6;
      b[8] = q5;
      b[9] = q4;
      b[10] = q3;
      b[11] = q2;
      a[3] = q7;
      b[12] = q6;
      b[13] = q5;
      b[14] = q4;
      b[15] = q3;
      int_math::gauss<4>(b, a, modulus);
      P.a[0] = a[0];
      P.a[1] = a[1];
      P.a[2] = a[2];
      P.a[3] = a[3];
      S.r[0] = q3;
      S.r[1] = q2;
      S.r[2] = q1;
      S.r[3] = q0;
      for (int i{0}; i < 4; ++i)
        backward();
    }
  }

  TRNG_CUDA_ENABLE
  inline void yarn4::jump2(unsigned int s) {
    result_type b[16], c[16]{};
    b[0] = P.a[0];
    b[1] = P.a[1];
    b[2] = P.a[2];
    b[3] = P.a[3];
    b[4] = 1;
    b[5] = 0;
    b[6] = 0;
    b[7] = 0;
    b[8] = 0;
    b[9] = 1;
    b[10] = 0;
    b[11] = 0;
    b[12] = 0;
    b[13] = 0;
    b[14] = 1;
    b[15] = 0;
    for (unsigned int i{0}; i < s; ++i) {
      int_math::matrix_mult<4>(b, b, c, modulus);
      ++i;
      if (not(i < s))
        break;
      int_math::matrix_mult<4>(c, c, b, modulus);
    }
    const result_type r[4]{S.r[0], S.r[1], S.r[2], S.r[3]};
    result_type d[4];
    if ((s & 1u) == 0)
      int_math::matrix_vec_mult<4>(b, r, d, modulus);
    else
      int_math::matrix_vec_mult<4>(c, r, d, modulus);
    S.r[0] = d[0];
    S.r[1] = d[1];
    S.r[2] = d[2];
    S.r[3] = d[3];
  }

  TRNG_CUDA_ENABLE
  inline void yarn4::jump(unsigned long long s) {
    if (s < 16) {
      for (unsigned int i{0}; i < s; ++i)
        step();
    } else {
      unsigned int i{0};
      while (s > 0) {
        if (s % 2 == 1)
          jump2(i);
        ++i;
        s >>= 1;
      }
    }
  }

  TRNG_CUDA_ENABLE
  inline void yarn4::discard(unsigned long long n) { jump(n); }

  TRNG_CUDA_ENABLE
  inline void yarn4::backward() {
    result_type t;
    if (P.a[3] != 0) {
      t = S.r[0];
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[0]) * static_cast<int64_t>(S.r[1])) % modulus);
      if (t < 0)
        t += modulus;
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[1]) * static_cast<int64_t>(S.r[2])) % modulus);
      if (t < 0)
        t += modulus;
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[2]) * static_cast<int64_t>(S.r[3])) % modulus);
      if (t < 0)
        t += modulus;
      t = static_cast<result_type>(
          (static_cast<int64_t>(t) *
           static_cast<int64_t>(int_math::modulo_invers(P.a[3], modulus))) %
          modulus);
    } else if (P.a[2] != 0) {
      t = S.r[1];
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[0]) * static_cast<int64_t>(S.r[2])) % modulus);
      if (t < 0)
        t += modulus;
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[1]) * static_cast<int64_t>(S.r[3])) % modulus);
      if (t < 0)
        t += modulus;
      t = static_cast<result_type>(
          (static_cast<int64_t>(t) *
           static_cast<int64_t>(int_math::modulo_invers(P.a[2], modulus))) %
          modulus);
    } else if (P.a[1] != 0) {
      t = S.r[2];
      t -= static_cast<result_type>(
          (static_cast<int64_t>(P.a[0]) * static_cast<int64_t>(S.r[3])) % modulus);
      if (t < 0)
        t += modulus;
      t = static_cast<result_type>(
          (static_cast<int64_t>(t) *
           static_cast<int64_t>(int_math::modulo_invers(P.a[1], modulus))) %
          modulus);
    } else if (P.a[0] != 0) {
      t = S.r[3];
      t = static_cast<result_type>(
          (static_cast<int64_t>(t) *
           static_cast<int64_t>(int_math::modulo_invers(P.a[0], modulus))) %
          modulus);
    } else
      t = 0;
    S.r[0] = S.r[1];
    S.r[1] = S.r[2];
    S.r[2] = S.r[3];
    S.r[3] = t;
  }

}  // namespace trng

#endif
