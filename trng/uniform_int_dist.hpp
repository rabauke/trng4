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

#if !(defined TRNG_UNIFORM_INT_DIST_HPP)

#define TRNG_UNIFORM_INT_DIST_HPP

#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <ostream>
#include <istream>
#include <ciso646>

namespace trng {

  // uniform random number generator class
  class uniform_int_dist {
  public:
    using result_type = int;

    class param_type {
    private:
      result_type a_{0}, b_{1}, d_{1};
      TRNG_CUDA_ENABLE
      result_type d() const { return d_; }

    public:
      TRNG_CUDA_ENABLE
      result_type a() const { return a_; }
      TRNG_CUDA_ENABLE
      void a(result_type a_new) {
        a_ = a_new;
        d_ = b_ - a_;
      }
      TRNG_CUDA_ENABLE
      result_type b() const { return b_; }
      TRNG_CUDA_ENABLE
      void b(result_type b_new) {
        b_ = b_new;
        d_ = b_ - a_;
      }
      TRNG_CUDA_ENABLE
      param_type() = default;
      TRNG_CUDA_ENABLE
      explicit param_type(result_type a, result_type b) : a_(a), b_(b), d_(b - a) {}

      friend class uniform_int_dist;
    };

  private:
    param_type P;

  public:
    // constructor
    TRNG_CUDA_ENABLE
    explicit uniform_int_dist(result_type a, result_type b) : P{a, b} {}
    TRNG_CUDA_ENABLE
    explicit uniform_int_dist(const param_type &P) : P{P} {}
    // reset internal state
    TRNG_CUDA_ENABLE
    void reset() {}
    // random numbers
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r) {
      return static_cast<result_type>(P.d() * utility::uniformco<double>(r)) + P.a();
    }
    template<typename R>
    TRNG_CUDA_ENABLE result_type operator()(R &r, const param_type &P) {
      uniform_int_dist g(P);
      return g(r);
    }
    // property methods
    TRNG_CUDA_ENABLE
    result_type min() const { return P.a(); }
    TRNG_CUDA_ENABLE
    result_type max() const { return P.b() - 1; }
    TRNG_CUDA_ENABLE
    param_type param() const { return P; }
    TRNG_CUDA_ENABLE
    void param(const param_type &p_new) { P = p_new; }
    TRNG_CUDA_ENABLE
    result_type a() const { return P.a(); }
    TRNG_CUDA_ENABLE
    void a(result_type a_new) { P.a(a_new); }
    TRNG_CUDA_ENABLE
    result_type b() const { return P.b(); }
    TRNG_CUDA_ENABLE
    void b(result_type b_new) { P.b(b_new); }
    // probability density function
    TRNG_CUDA_ENABLE
    double pdf(result_type x) const {
      if (x < P.a() or x >= P.b())
        return 0.0;
      return 1.0 / P.d();
    }
    // cumulative density function
    TRNG_CUDA_ENABLE
    double cdf(result_type x) const {
      if (x < P.a())
        return 0;
      if (x >= P.b())
        return 1.0;
      return static_cast<double>(x - P.a() + 1) / static_cast<double>(P.d());
    }
  };

  // -------------------------------------------------------------------

  // EqualityComparable concept
  TRNG_CUDA_ENABLE
  bool operator==(const uniform_int_dist::param_type &P1,
                  const uniform_int_dist::param_type &P2) {
    return P1.a() == P2.a() and P1.b() == P2.b();
  }
  TRNG_CUDA_ENABLE
  bool operator!=(const uniform_int_dist::param_type &P1,
                  const uniform_int_dist::param_type &P2) {
    return not(P1 == P2);
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const uniform_int_dist::param_type &P) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << '(' << P.a() << ' ' << P.b() << ')';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   uniform_int_dist::param_type &P) {
    int a, b;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::delim('(') >> a >> utility::delim(' ') >> b >> utility::delim(')');
    if (in)
      P = uniform_int_dist::param_type(a, b);
    in.flags(flags);
    return in;
  }

  // -------------------------------------------------------------------

  // EqualityComparable concept
  TRNG_CUDA_ENABLE
  bool operator==(const uniform_int_dist &g1, const uniform_int_dist &g2) {
    return g1.param() == g2.param();
  }
  TRNG_CUDA_ENABLE
  bool operator!=(const uniform_int_dist &g1, const uniform_int_dist &g2) {
    return g1.param() != g2.param();
  }

  // Streamable concept
  template<typename char_t, typename traits_t>
  std::basic_ostream<char_t, traits_t> &operator<<(std::basic_ostream<char_t, traits_t> &out,
                                                   const uniform_int_dist &g) {
    std::ios_base::fmtflags flags(out.flags());
    out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    out << "[uniform_int " << g.param() << ']';
    out.flags(flags);
    return out;
  }

  template<typename char_t, typename traits_t>
  std::basic_istream<char_t, traits_t> &operator>>(std::basic_istream<char_t, traits_t> &in,
                                                   uniform_int_dist &g) {
    uniform_int_dist::param_type P;
    std::ios_base::fmtflags flags(in.flags());
    in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
    in >> utility::ignore_spaces() >> utility::delim("[uniform_int ") >> P >>
        utility::delim(']');
    if (in)
      g.param(P);
    in.flags(flags);
    return in;
  }

}  // namespace trng

#endif
