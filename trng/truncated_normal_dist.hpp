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

#if !(defined TRNG_TRUNCATED_NORMAL_DIST_HPP)

#define TRNG_TRUNCATED_NORMAL_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/constants.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <ciso646>

namespace trng {

  // uniform random number generator class
  template<typename float_t = double>
  class truncated_normal_dist {
  public:
    typedef float_t result_type;
    class param_type;

    class param_type {
    private:
      result_type mu_, sigma_, a_, b_, Phi_a, Phi_b;

      TRNG_CUDA_ENABLE
      void update_Phi() {
        if (a_ != -math::numeric_limits<result_type>::infinity())
          Phi_a = math::Phi((a_ - mu_) / sigma_);
        else
          Phi_a = result_type(0);
        if (b_ != math::numeric_limits<result_type>::infinity())
          Phi_b = math::Phi((b_ - mu_) / sigma_);
        else
          Phi_b = result_type(1);
      }

    public:
      TRNG_CUDA_ENABLE
      result_type mu() const { return mu_; }
      TRNG_CUDA_ENABLE
      void mu(result_type mu_new) {
        mu_ = mu_new;
        update_Phi();
      }
      TRNG_CUDA_ENABLE
      result_type sigma() const { return sigma_; }
      TRNG_CUDA_ENABLE
      void sigma(result_type sigma_new) {
        sigma_ = sigma_new;
        update_Phi();
      }
      TRNG_CUDA_ENABLE
      result_type a() const { return a_; }
      TRNG_CUDA_ENABLE
      void a(result_type a_new) {
        a_ = a_new;
        update_Phi();
      }
      TRNG_CUDA_ENABLE
      result_type b() const { return b_; }
      TRNG_CUDA_ENABLE
      void b(result_type b_new) {
        b_ = b_new;
        update_Phi();
      }
      TRNG_CUDA_ENABLE
      param_type()
          : mu_(0),
            sigma_(1),
            a_(-math::numeric_limits<result_type>::infinity()),
            b_(math::numeric_limits<result_type>::infinity()) {
        update_Phi();
      }
      TRNG_CUDA_ENABLE
      param_type(result_type mu, result_type sigma, result_type a, result_type b)
          : mu_(mu), sigma_(sigma), a_(a), b_(b) {
        update_Phi();
      }

      friend class truncated_normal_dist;

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const param_type &p) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << std::setprecision(math::numeric_limits<float_t>::digits10 + 1) << p.mu()
            << ' ' << p.sigma() << ' ' << p.a() << ' ' << p.b() << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, param_type &p) {
        float_t mu, sigma, a, b;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> mu >> utility::delim(' ') >> sigma >>
            utility::delim(' ') >> a >> utility::delim(' ') >> b >> utility::delim(')');
        if (in)
          p = param_type(mu, sigma, a, b);
        in.flags(flags);
        return in;
      }
    };

  private:
    param_type p;

  public:
    // constructor
    TRNG_CUDA_ENABLE
    truncated_normal_dist(result_type mu, result_type sigma, result_type a, result_type b)
        : p(mu, sigma, a, b) {}
    TRNG_CUDA_ENABLE
    explicit truncated_normal_dist(const param_type &p) : p(p) {}
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() {}
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r) {
      return icdf(utility::uniformoo<result_type>(r));
    }
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r, const param_type &p) {
      truncated_normal_dist g(p);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return p.a(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return p.b(); }
    TRNG_CUDA_ENABLE
    param_type param() const { return p; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { p = p_new; }
    TRNG_CUDA_ENABLE
    result_type mu() const { return p.mu(); }
    TRNG_CUDA_ENABLE
    void mu(result_type mu_new) { p.mu(mu_new); }
    TRNG_CUDA_ENABLE
    result_type sigma() const { return p.sigma(); }
    TRNG_CUDA_ENABLE
    void sigma(result_type sigma_new) { p.sigma(sigma_new); }
    TRNG_CUDA_ENABLE
    result_type a() const { return p.a(); }
    TRNG_CUDA_ENABLE
    void a(result_type a_new) { p.a(a_new); }
    TRNG_CUDA_ENABLE
    result_type b() const { return p.b(); }
    TRNG_CUDA_ENABLE
    void b(result_type b_new) { p.b(b_new); }
    // probability density function
    TRNG_CUDA_ENABLE
    result_type pdf(result_type x) const {
      x -= p.mu();
      x /= p.sigma();
      return math::constants<result_type>::one_over_sqrt_2pi() / p.sigma() *
             math::exp(-0.5 * x * x) / (p.Phi_b - p.Phi_a);
    }
    // cumulative density function
    TRNG_CUDA_ENABLE
    result_type cdf(result_type x) const {
      x -= p.mu();
      x /= p.sigma();
      return (math::Phi(x) - p.Phi_a) / (p.Phi_b - p.Phi_a);
    }
    // inverse cumulative density function
    TRNG_CUDA_ENABLE
    result_type icdf(result_type x) const {
      x *= p.Phi_b - p.Phi_a;
      x += p.Phi_a;
      return math::inv_Phi(x) * p.sigma() + p.mu();
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(
      const typename truncated_normal_dist<float_t>::param_type &p1,
      const typename truncated_normal_dist<float_t>::param_type &p2) {
    return p1.mu() == p2.mu() and p1.sigma() == p2.sigma();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(
      const typename truncated_normal_dist<float_t>::param_type &p1,
      const typename truncated_normal_dist<float_t>::param_type &p2) {
    return not(p1 == p2);
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator==(const truncated_normal_dist<float_t> &g1,
                                          const truncated_normal_dist<float_t> &g2) {
    return g1.param() == g2.param();
  }

  template<typename float_t>
  TRNG_CUDA_ENABLE inline bool operator!=(const truncated_normal_dist<float_t> &g1,
                                          const truncated_normal_dist<float_t> &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const truncated_normal_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[truncated_normal " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   truncated_normal_dist<float_t> &g) {
    typename truncated_normal_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[truncated_normal ") >> p >>
        utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
