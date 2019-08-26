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

#if !(defined TRNG_CORRELATED_NORMAL_DIST_HPP)

#define TRNG_CORRELATED_NORMAL_DIST_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/math.hpp>
#include <trng/special_functions.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <cerrno>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ciso646>

namespace trng {

  // uniform random number generator class
  template<typename float_t = double>
  class correlated_normal_dist {
  public:
    typedef float_t result_type;
    class param_type;

    class param_type {
      typedef typename std::vector<result_type>::size_type size_type;
      std::vector<result_type> H_;
      size_type d_{0};

      void Cholesky_factorization() {
        for (size_type i = 0; i < d_; ++i) {
          result_type t;
          for (size_type k = 0; k < i; ++k) {
            t = 0.0;
            for (size_type j = 0; j < k; ++j)
              t += H_[i * d_ + j] * H_[k * d_ + j];
            H_[i * d_ + k] = (H_[i * d_ + k] - t) / H_[k * d_ + k];
          }
          t = 0.0;
          for (size_type j = 0; j < i; ++j)
            t += H_[i * d_ + j] * H_[i * d_ + j];
          H_[i * d_ + i] = trng::math::sqrt(H_[i * d_ + i] - t);
          for (size_type k = i + 1; k < d_; ++k)
            H_[i * d_ + k] = 0;
        }
      }

      result_type H_times(const std::vector<result_type> &normal) {
        size_type d(normal.size());
        result_type y = 0.0;
        for (size_type j = 0; j < d; ++j)
          y += H_[(d - 1) * d_ + j] * normal[j];
        return y;
      }

      size_type d() const { return d_; }

    public:
      param_type() = default;
      template<typename iter>
      explicit param_type(iter first, iter last) {
        d_ = static_cast<size_type>(
            trng::math::sqrt(static_cast<result_type>(std::distance(first, last))));
        H_.reserve(d_ * d_);
        for (size_type i = 0; i < d_ * d_; ++i, ++first)
          H_.push_back(*first);
        Cholesky_factorization();
      }

      friend class correlated_normal_dist;

      // EqualityComparable concept
      friend inline bool operator==(const typename correlated_normal_dist::param_type &p1,
                                    const typename correlated_normal_dist::param_type &p2) {
        return p1.H_ == p2.H_;
      }
      friend inline bool operator!=(const typename correlated_normal_dist::param_type &p1,
                                    const typename correlated_normal_dist::param_type &p2) {
        return not(p1 == p2);
      }

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const param_type &p) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << p.d_ << std::setprecision(17);
        for (unsigned int i = 0; i < p.d_ * p.d_; ++i)
          out << ' ' << p.H_[i];
        out << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, param_type &p) {
        std::vector<result_type> H;
        unsigned int d;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> d;
        for (unsigned int i = 0; i < d * d; ++i) {
          result_type t = 0;
          in >> utility::delim(' ') >> t;
          H.push_back(t);
        }
        in >> utility::delim(')');
        if (in) {
          p.d_ = d;
          p.H_ = H;
        }
        in.flags(flags);
        return in;
      }
    };

  private:
    param_type p;
    std::vector<result_type> normal;

  public:
    // constructor
    template<typename iter>
    correlated_normal_dist(iter first, iter last) : p(first, last) {}
    explicit correlated_normal_dist(const param_type &p) : p(p) {}
    // reset internal state
    void reset() { normal.clear(); }
    // random numbers
    template<typename R>
    result_type operator()(R &r) {
      normal.push_back(trng::math::inv_Phi(utility::uniformoo<result_type>(r)));
      result_type y = p.H_times(normal);
      if (normal.size() == p.d())
        normal.clear();
      return y;
    }
    template<typename R>
    result_type operator()(R &r, const param_type &p) {
      correlated_normal_dist g(p);
      return g(r);
    }
    // property methods
    result_type min() const { return -math::numeric_limits<result_type>::infinity(); }
    result_type max() const { return math::numeric_limits<result_type>::infinity(); }
    param_type param() const { return p; }
    void param(const param_type &p_new) { p = p_new; }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  template<typename float_t>
  inline bool operator==(const correlated_normal_dist<float_t> &g1,
                         const correlated_normal_dist<float_t> &g2) {
    return g1.param() == g2.param();
  }
  template<typename float_t>
  inline bool operator!=(const correlated_normal_dist<float_t> &g1,
                         const correlated_normal_dist<float_t> &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t, typename float_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const correlated_normal_dist<float_t> &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[correlated_normal " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t, typename float_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   correlated_normal_dist<float_t> &g) {
    typename correlated_normal_dist<float_t>::param_type p;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[correlated_normal ") >> p >>
        utility::delim(']');
    if (in)
      g.param(p);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
