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

#if !(defined TRNG_STUDENT_T_DIST_HPP)

#define TRNG_STUDENT_T_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>
#include <ciso646>

namespace trng {

  // uniform random number generator class
  template<typename float_t = double>
  class student_t_dist {
  public:
    using result_type = float_t;

    class param_type {
    private:
      int nu_{1};

    public:
      TRNG_CUDA_ENABLE
      int nu() const { return nu_; }
      TRNG_CUDA_ENABLE
      void nu(int nu_new) { nu_ = nu_new; }
      TRNG_CUDA_ENABLE
      param_type() = default;
      TRNG_CUDA_ENABLE
      explicit param_type(int nu) : nu_(nu) {}

      friend class student_t_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const param_type &P) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << P.nu() << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, param_type &P) {
        int nu;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> nu >> utility::delim(')');
        if (in)
          P = param_type(nu);
        in.flags(flags);
        return in;
      }
    };

  private:
    param_type P;

    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      const result_type t{
          math::inv_Beta_I(x, P.nu() / result_type(2), P.nu() / result_type(2))};
      return math::sqrt(P.nu() / (t * (1 - t))) * (t - result_type(1) / result_type(2));
    }

  public:
    // constructor
    TRNG_CUDA_ENABLE
    explicit student_t_dist(int nu) : P{nu} {}
    TRNG_CUDA_ENABLE
    explicit student_t_dist(const param_type &P) : P{P} {}
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() {}
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r) {
      return icdf_(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r, const param_type &P) {
      student_t_dist g(P);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return -math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return P; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { P = p_new; }
    TRNG_CUDA_ENABLE
    int nu() const { return P.nu(); }
    TRNG_CUDA_ENABLE
    void nu(int nu_new) { P.nu(nu_new); }
    // probability density function
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      const result_type norm{math::exp(math::ln_Gamma((P.nu() + 1) / result_type(2)) -
                                       math::ln_Gamma(P.nu() / result_type(2))) /
                             math::sqrt(math::constants<result_type>::pi * P.nu())};
      return norm * math::pow(1 + x * x / P.nu(), (P.nu() + 1) / result_type(-2));
    }
    // cumulative density function
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      const result_type t1{+math::sqrt(x * x + P.nu())};
      const result_type t2{(x + t1) / (2 * t1)};
      return math::Beta_I(t2, P.nu() / result_type(2), P.nu() / result_type(2));
    }
    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf(result_type x) const {
      if (x <= 0 or x >= 1) {
#if !(defined __CUDA_ARCH__)
        errno = EDOM;
#endif
        return math::numeric_limits<result_type>::quiet_NaN();
      }
      if (x == 0)
        return -math::numeric_limits<result_type>::infinity();
      if (x == 1)
        return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(
      const typename student_t_dist<float_t>::param_type &P1,
      const typename student_t_dist<float_t>::param_type &P2) {
    return P1.nu() == P2.nu();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(
      const typename student_t_dist<float_t>::param_type &P1,
      const typename student_t_dist<float_t>::param_type &P2) {
    return not(P1 == P2);
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(const student_t_dist<float_t> &g1,
                                          const student_t_dist<float_t> &g2) {
    return g1.param() == g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(const student_t_dist<float_t> &g1,
                                          const student_t_dist<float_t> &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const student_t_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[student_t " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   student_t_dist<float_t> &g) {
    typename student_t_dist<float_t>::param_type P;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[student_t ") >> P >> utility::delim(']');
    if (in)
      g.param(P);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
