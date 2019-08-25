// Copyright (c) 2000-2019, Heiko Bauke
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

#if !(defined TRNG_SNEDECOR_F_DIST_HPP)

#define TRNG_SNEDECOR_F_DIST_HPP

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
  class snedecor_f_dist {
  public:
    typedef float_t result_type;
    class param_type;

    class param_type {
    private:
      int n_, m_;

    public:
      TRNG_CUDA_ENABLE
      int n() const { return n_; }
      TRNG_CUDA_ENABLE
      void n(int n_new) { n_ = n_new; }
      TRNG_CUDA_ENABLE
      int m() const { return m_; }
      TRNG_CUDA_ENABLE
      void m(int m_new) { m_ = m_new; }
      TRNG_CUDA_ENABLE
      param_type() : n_(1), m_(1) {}
      TRNG_CUDA_ENABLE
      param_type(int n, int m) : n_(n), m_(m) {}

      friend class snedecor_f_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const param_type &p) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << p.n() << ' ' << p.m() << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, param_type &p) {
        int n, m;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> n >> utility::delim(' ') >> m >> utility::delim(')');
        if (in)
          p = snedecor_f_dist::param_type(n, m);
        in.flags(flags);
        return in;
      }
    };

  private:
    param_type p;

    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf_(result_type x) const {
      result_type t = math::inv_Beta_I(x, result_type(1) / result_type(2) * p.n(),
                                       result_type(1) / result_type(2) * p.m());
      return t / (1 - t) * static_cast<result_type>(p.m()) / static_cast<result_type>(p.n());
    }

  public:
    // constructor
    TRNG_CUDA_ENABLE
    snedecor_f_dist(int n, int m) : p(n, m) {}
    TRNG_CUDA_ENABLE
    explicit snedecor_f_dist(const param_type &p) : p(p) {}
    // reset internal state
    void reset() {}
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r) {
      return icdf_(utility::uniformco<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r, const param_type &p) {
      snedecor_f_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return 0; }
    TRNG_CUDA_ENABLE
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p = p_new; }
    TRNG_CUDA_ENABLE
    int n() const { return p.n(); }
    TRNG_CUDA_ENABLE
    void n(int n_new) { p.n(n_new); }
    TRNG_CUDA_ENABLE
    int m() const { return p.m(); }
    TRNG_CUDA_ENABLE
    void m(int m_new) { p.m(m_new); }
    // probability density function
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      const result_type n = static_cast<result_type>(p.n()),
                        m = static_cast<result_type>(p.m());
      // return math::exp(math::ln(n)*result_type(1)/result_type(2)*n +
      // math::ln(m)*result_type(1)/result_type(2)*m -
      // 		       math::ln(m+n*x)*result_type(1)/result_type(2)*(m+n) +
      // 		       math::ln(x)*(result_type(1)/result_type(2)*n-1))/math::Beta(result_type(1)/result_type(2)*n,
      // result_type(1)/result_type(2)*m);
      return math::exp(result_type(1) / result_type(2) * math::ln(n / m) * n +
                       math::ln(x) * (result_type(1) / result_type(2) * n - 1) -
                       math::ln(1 + n * x / m) * (result_type(1) / result_type(2) * n +
                                                  result_type(1) / result_type(2) * m) -
                       math::ln_Gamma(result_type(1) / result_type(2) * n) -
                       math::ln_Gamma(result_type(1) / result_type(2) * m) +
                       math::ln_Gamma(result_type(1) / result_type(2) * (n + m)));
    }
    // cumulative density function
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      const result_type n = static_cast<result_type>(p.n()),
                        m = static_cast<result_type>(p.m());
      return math::Beta_I(n * x / (m + n * x), result_type(1) / result_type(2) * n,
                          result_type(1) / result_type(2) * m);
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
        return 0;
      if (x == 1)
        return math::numeric_limits<result_type>::infinity();
      return icdf_(x);
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(
      const typename snedecor_f_dist<float_t>::param_type &p1,
      const typename snedecor_f_dist<float_t>::param_type &p2) {
    return p1.n() == p2.n() and p1.m() == p2.m();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(
      const typename snedecor_f_dist<float_t>::param_type &p1,
      const typename snedecor_f_dist<float_t>::param_type &p2) {
    return not(p1 == p2);
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(const snedecor_f_dist<float_t> &g1,
                                          const snedecor_f_dist<float_t> &g2) {
    return g1.param() == g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(const snedecor_f_dist<float_t> &g1,
                                          const snedecor_f_dist<float_t> &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const snedecor_f_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[snedecor_f " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   snedecor_f_dist<float_t> &g) {
    typename snedecor_f_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[snedecor_f ") >> p >>
        utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
