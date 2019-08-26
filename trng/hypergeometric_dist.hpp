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

#if !(defined TRNG_HYPERGEOMETRIC_DIST_HPP)

#define TRNG_HYPERGEOMETRIC_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <ciso646>

namespace trng {

  // non-uniform random number generator class
  class hypergeometric_dist {
  public:
    typedef int result_type;
    class param_type;

    class param_type {
    private:
      int n_{0},               // total number of balls in urn
          m_{0},               // number of "white" balls in urn
          d_{0},               // number of selected balls
          x_min{0}, x_max{0};  // minimum and maximum values of random variable
      std::vector<double> P_;

      void calc_probabilities() {
        x_min = std::max(0, d_ - n_ + m_);
        x_max = std::min(d_, m_);
        P_ = std::vector<double>();
        for (int x = x_min; x <= x_max; ++x)
          P_.push_back(math::exp(
              math::ln_binomial(static_cast<double>(m_), static_cast<double>(x)) +
              math::ln_binomial(static_cast<double>(n_ - m_), static_cast<double>(m_ - x))));
        // build list with cumulative density function
        for (std::vector<double>::size_type i(1); i < P_.size(); ++i)
          P_[i] += P_[i - 1];
        for (std::vector<double>::size_type i(0); i < P_.size(); ++i)
          P_[i] /= P_.back();
      }

    public:
      int n() const { return n_; }
      void n(int n_new) {
        n_ = n_new;
        calc_probabilities();
      }
      int m() const { return m_; }
      void m(int m_new) {
        m_ = m_new;
        calc_probabilities();
      }
      int d() const { return d_; }
      void d(int d_new) {
        d_ = d_new;
        calc_probabilities();
      }
      param_type() = default;
      explicit param_type(int n, int m, int d) : n_(n), m_(m), d_(d) { calc_probabilities(); }
      friend class hypergeometric_dist;
    };

  private:
    param_type P;

  public:
    // constructor
    hypergeometric_dist(int n, int m, int d) : P(n, m, d) {}
    explicit hypergeometric_dist(const param_type &P) : P(P) {}
    // reset internal state
    void reset() {}
    // random numbers
    template<typename R>
    int operator()(R &r) {
      return P.x_min +
             utility::discrete(utility::uniformoo<double>(r), P.P_.begin(), P.P_.end());
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      hypergeometric_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return P.x_min; }
    int max() const { return P.x_max; }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P = P_new; }
    int n() const { return P.n(); }
    void n(int n_new) { P.n(n_new); }
    int m() const { return P.m(); }
    void m(int m_new) { P.m(m_new); }
    int d() const { return P.d(); }
    void d(int d_new) { P.d(d_new); }
    // probability density function
    double pdf(int x) const {
      if (x < P.x_min or x > P.x_max)
        return 0.0;
      x -= P.x_min;
      if (x == 0)
        return P.P_[0];
      return P.P_[x] - P.P_[x - 1];
    }
    // cumulative density function
    double cdf(int x) const {
      if (x < P.x_min)
        return 0.0;
      if (x > P.x_max)
        return 1.0;
      return P.P_[x - P.x_min];
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const hypergeometric_dist::param_type &p1,
                         const hypergeometric_dist::param_type &p2) {
    return p1.n() == p2.n() and p1.m() == p2.m() and p1.d() == p2.d();
  }
  inline bool operator!=(const hypergeometric_dist::param_type &p1,
                         const hypergeometric_dist::param_type &p2) {
    return not(p1 == p2);
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const hypergeometric_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << '(' << std::setprecision(17) << P.n() << ' ' << P.m() << ' ' << P.d() << ')';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   hypergeometric_dist::param_type &P) {
    int n, m, d;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::delim('(') >> n >> utility::delim(' ') >> m >> utility::delim(' ') >> d >>
        utility::delim(')');
    if (in)
      P = hypergeometric_dist::param_type(n, m, d);
    in.flags(flags);
    return in;
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const hypergeometric_dist &g1, const hypergeometric_dist &g2) {
    return g1.param() == g2.param();
  }
  inline bool operator!=(const hypergeometric_dist &g1, const hypergeometric_dist &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const hypergeometric_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[hypergeometric " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   hypergeometric_dist &g) {
    hypergeometric_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[hypergeometric ") >> p >>
        utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
