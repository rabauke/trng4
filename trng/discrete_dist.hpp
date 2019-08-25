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
//     with the disctribution.
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

#if !(defined TRNG_DISCRETE_DIST_HPP)

#define TRNG_DISCRETE_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <ostream>
#include <iomanip>
#include <istream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ciso646>

namespace trng {

  // non-uniform random number generator class
  class discrete_dist {
  public:
    typedef int result_type;
    class param_type;

    class param_type {
    private:
      typedef std::vector<double>::size_type size_type;
      std::vector<double> P;
      size_type N, offset, layers;

    public:
      param_type() : P(), N(0), offset(0), layers(0) {}
      template<typename iter>
      param_type(iter first, iter last) : P(first, last) {
        N = P.size();
        layers = math::log2_ceil(N);
        offset = math::pow2(layers) - 1;
        P.resize(N + offset);
        std::copy_backward(P.begin(), P.begin() + N, P.end());
        std::fill(P.begin(), P.begin() + offset, 0);
        update_all_layers();
      }
      explicit param_type(int n) : P(n, 1.0) {
        N = P.size();
        layers = math::log2_ceil(N);
        offset = math::pow2(layers) - 1;
        P.resize(N + offset);
        std::copy_backward(P.begin(), P.begin() + N, P.end());
        std::fill(P.begin(), P.begin() + offset, 0);
        update_all_layers();
      }

    private:
      void update_layer(size_type layer, size_type n) {
        size_type first = math::pow2(layer) - 1, last = first + n;
        for (size_type i = first; i < last; ++i, ++i)
          if (i + 1 < last)
            P[(i - 1) / 2] = P[i] + P[i + 1];
          else
            P[(i - 1) / 2] = P[i];
      }
      void update_all_layers() {
        size_type layer = layers;
        if (layer > 0) {
          update_layer(layer, N);
          --layer;
        }
        while (layer > 0) {
          update_layer(layer, math::pow2(layer));
          --layer;
        }
      }

    public:
      friend class discrete_dist;
      friend bool operator==(const param_type &, const param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &, const discrete_dist::param_type &);
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &, discrete_dist::param_type &);
    };

  private:
    param_type P;

  public:
    // constructor
    template<typename iter>
    discrete_dist(iter first, iter last) : P(first, last) {}
    explicit discrete_dist(int N) : P(N) {}
    explicit discrete_dist(const param_type &P) : P(P) {}
    // reset internal state
    void reset() {}
    // random numbers
    template<typename R>
    int operator()(R &r) {
      if (P.N == 0)
        return -1;
      double u(utility::uniformco<double>(r) * P.P[0]);
      param_type::size_type x(0);
      while (x < P.offset) {
        if (u < P.P[2 * x + 1]) {
          x = 2 * x + 1;
        } else {
          u -= P.P[2 * x + 1];
          x = 2 * x + 2;
        }
      }
      return static_cast<int>(x - P.offset);
    }
    template<typename R>
    int operator()(R &r, const param_type &p) {
      discrete_dist g(p);
      return g(r);
    }
    // property methods
    int min() const { return 0; }
    int max() const { return static_cast<int>(P.N - 1); }
    param_type param() const { return P; }
    void param(const param_type &P_new) { P = P_new; }
    void param(int x, double p) {
      x += static_cast<int>(P.offset);
      P.P[x] = p;
      if (x > 0) {
        do {
          x = (x - 1) / 2;
          P.P[x] = P.P[2 * x + 1] + P.P[2 * x + 2];
        } while (x > 0);
      }
    }
    // probability density function
    double pdf(int x) const {
      return (x < 0 or x >= static_cast<int>(P.N)) ? 0.0 : P.P[x + P.offset] / P.P[0];
    }
    // cumulative density function
    double cdf(int x) const {
      if (x < 0)
        return 0.0;
      if (x < static_cast<int>(P.N))
        return std::accumulate(&P.P[P.offset], &P.P[x + P.offset + 1], 0.0) / P.P[0];
      return 1.0;
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const discrete_dist::param_type &p1,
                         const discrete_dist::param_type &p2) {
    return p1.P == p2.P;
  }
  inline bool operator!=(const discrete_dist::param_type &p1,
                         const discrete_dist::param_type &p2) {
    return !(p1 == p2);
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const discrete_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << '(' << P.N << ' ';
    for (std::vector<double>::size_type i = P.offset; i < P.P.size(); ++i) {
      out << std::setprecision(17) << P.P[i];
      if (i + 1 < P.P.size())
        out << ' ';
    }
    out << ')';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   discrete_dist::param_type &P) {
    double p;
    std::vector<double>::size_type n;
    std::vector<double> P_new;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::delim('(') >> n >> utility::delim(' ');
    for (std::vector<double>::size_type i = 0; i < n; ++i) {
      in >> p;
      if (i + 1 < n)
        in >> utility::delim(' ');
      P_new.push_back(p);
    }
    in >> utility::delim(')');
    if (in)
      P = discrete_dist::param_type(P_new.begin(), P_new.end());
    in.flags(flags);
    return in;
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  inline bool operator==(const discrete_dist &g1, const discrete_dist &g2) {
    return g1.param() == g2.param();
  }
  inline bool operator!=(const discrete_dist &g1, const discrete_dist &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const discrete_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[discrete " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   discrete_dist &g) {
    discrete_dist::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[discrete ") >> p >> utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
