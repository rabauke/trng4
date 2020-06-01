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

#if !(defined TRNG_CAUCHY_DIST_HPP)

#define TRNG_CAUCHY_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/constants.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>
#include <ciso646>

namespace trng {

  // uniform random number generator class
  template<typename float_t = double>
  class cauchy_dist {
  public:
    using result_type = float_t;

    class param_type {
    private:
      result_type theta_{1}, eta_{0};

    public:
      TRNG_CUDA_ENABLE
      result_type theta() const { return theta_; }
      TRNG_CUDA_ENABLE
      void theta(result_type theta_new) { theta_ = theta_new; }
      TRNG_CUDA_ENABLE
      result_type eta() const { return eta_; }
      TRNG_CUDA_ENABLE
      void eta(result_type eta_new) { eta_ = eta_new; }
      TRNG_CUDA_ENABLE
      param_type() = default;
      TRNG_CUDA_ENABLE
      explicit param_type(result_type theta, result_type eta) : theta_(theta), eta_(eta) {}

      friend class cauchy_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const param_type &P) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << std::setprecision(math::numeric_limits<float_t>::digits10 + 1)
            << P.theta() << ' ' << P.eta() << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, param_type &P) {
        float_t theta, eta;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> theta >> utility::delim(' ') >> eta >> utility::delim(')');
        if (in)
          P = param_type(theta, eta);
        in.flags(flags);
        return in;
      }
    };

  private:
    param_type P;

    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      const result_type t{x * math::constants<result_type>::pi};
      return P.eta() - P.theta() * math::cos(t) / math::sin(t);
    }

  public:
    // constructor
    TRNG_CUDA_ENABLE
    explicit cauchy_dist(result_type theta, result_type eta) : P{theta, eta} {}
    TRNG_CUDA_ENABLE
    explicit cauchy_dist(const param_type &P) : P{P} {}
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
      cauchy_dist g(P);
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
    result_type theta() const { return P.theta(); }
    TRNG_CUDA_ENABLE
    void theta(result_type theta_new) { P.theta(theta_new); }
    TRNG_CUDA_ENABLE
    result_type eta() const { return P.eta(); }
    TRNG_CUDA_ENABLE
    void eta(result_type eta_new) { P.eta(eta_new); }
    // probability density function
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      x -= P.eta();
      x /= P.theta();
      return math::constants<result_type>::one_over_pi / (1 + x * x) / P.theta();
    }
    // cumulative density function
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      x -= P.eta();
      x /= P.theta();
      return math::constants<result_type>::one_over_pi * math::atan(x) +
             result_type(1) / result_type(2);
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
  TRNG_CUDA_ENABLE inline bool operator==(const typename cauchy_dist<float_t>::param_type &P1,
                                          const typename cauchy_dist<float_t>::param_type &P2) {
    return P1.theta() == P2.theta() and P1.eta() == P2.eta();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(const typename cauchy_dist<float_t>::param_type &P1,
                                          const typename cauchy_dist<float_t>::param_type &P2) {
    return not(P1 == P2);
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(const cauchy_dist<float_t> &g1,
                                          const cauchy_dist<float_t> &g2) {
    return g1.param() == g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(const cauchy_dist<float_t> &g1,
                                          const cauchy_dist<float_t> &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const cauchy_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[cauchy " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   cauchy_dist<float_t> &g) {
    typename cauchy_dist<float_t>::param_type P;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[cauchy ") >> P >> utility::delim(']');
    if (in)
      g.param(P);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
