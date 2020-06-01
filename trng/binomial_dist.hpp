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

#if !(defined TRNG_BINOMIAL_DIST_HPP)

#define TRNG_BINOMIAL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <vector>
#include <ciso646>

namespace trng {

  // non-uniform random number generator class
  class binomial_dist {
  public:
    using result_type = int;

    class param_type {
    private:
      double p_{0.5};
      int n_{0};
      std::vector<double> P_;

      void calc_probabilities() {
        P_ = std::vector<double>();
        P_.reserve(n_ + 1);
        double ln_binom{0.0};
        const double ln_p{math::ln(p_)};
        const double ln_1_p{math::ln(1.0 - p_)};
        for (int i{0}; i <= n_; ++i) {
          const double ln_prob{ln_binom + static_cast<double>(i) * ln_p +
                               static_cast<double>(n_ - i) * ln_1_p};
          P_.push_back(math::exp(ln_prob));
          ln_binom += math::ln(static_cast<double>(n_ - i));
          ln_binom -= math::ln(static_cast<double>(i + 1));
        }
        // build list with cumulative density function
        for (std::vector<double>::size_type i{1}; i < P_.size(); ++i)
          P_[i] += P_[i - 1];
        // normailze, just in case of rounding errors
        for (std::vector<double>::size_type i{0}; i < P_.size(); ++i)
          P_[i] /= P_.back();
      }

    public:
      double p() const { return p_; }
      void p(double p_new) {
        p_ = p_new;
        calc_probabilities();
      }
      int n() const { return n_; }
      void n(int n_new) {
        n_ = n_new;
        calc_probabilities();
      }
      param_type() = default;
      explicit param_type(double p, int n) : p_(p), n_(n) { calc_probabilities(); }
      friend class binomial_dist;
    };

  private:
    param_type P;

  public:
    // constructor
    explicit binomial_dist(double p, int n) : P{p, n} {}
    explicit binomial_dist(const param_type &P) : P{P} {}
    // reset internal state
    void reset() {}
    // random numbers
    template<typename R>
    int operator()(R &r) {
      return utility::discrete(utility::uniformoo<double>(r), P.P_.begin(), P.P_.end());
    }
    template<typename R>
    int operator()(R &r, const param_type &P) {
      binomial_dist g(P);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return P.n(); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P = P_new; }
    double p() const { return P.p(); }
    void p(double p_new) { P.p(p_new); }
    int n() const { return P.n(); }
    void n(int n_new) { P.n(n_new); }
    // probability density function
    double pdf(int x) const {
      if (x < 0 or x > P.n())
        return 0.0;
      if (x == 0)
        return P.P_[0];
      return P.P_[x] - P.P_[x - 1];
    }
    // cumulative density function
    double cdf(int x) const {
      if (x < 0)
        return 0.0;
      if (x <= P.n())
        return P.P_[x];
      return 1.0;
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const binomial_dist::param_type &P1,
                         const binomial_dist::param_type &P2) {
    return P1.p() == P2.p() and P1.n() == P2.n();
  }
  inline bool operator!=(const binomial_dist::param_type &P1,
                         const binomial_dist::param_type &P2) {
    return not(P1 == P2);
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const binomial_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << '(' << std::setprecision(17) << P.p() << ' ' << P.n() << ')';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   binomial_dist::param_type &P) {
    double p;
    int n;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::delim('(') >> p >> utility::delim(' ') >> n >> utility::delim(')');
    if (in)
      P = binomial_dist::param_type(p, n);
    in.flags(flags);
    return in;
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const binomial_dist &g1, const binomial_dist &g2) {
    return g1.param() == g2.param();
  }
  inline bool operator!=(const binomial_dist &g1, const binomial_dist &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const binomial_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[binomial " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   binomial_dist &g) {
    binomial_dist::param_type P;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[binomial ") >> P >> utility::delim(']');
    if (in)
      g.param(P);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
