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

#if !(defined TRNG_LAGFIB2XOR_HPP)

#define TRNG_LAGFIB2XOR_HPP

#include <trng/limits.hpp>
#include <trng/utility.hpp>
#include <trng/minstd.hpp>
#include <trng/int_types.hpp>
#include <trng/linear_algebra.hpp>
#include <climits>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <sstream>
#include <ciso646>

namespace trng {

  template<typename integer_type, unsigned int A, unsigned int B>
  class lagfib2xor {
  public:
    // Uniform random number generator concept
    using result_type = integer_type;
    result_type operator()() {
      step();
      return S.r[S.index];
    }

  private:
    static constexpr result_type min_ = 0;
    static constexpr result_type max_ = ~result_type(0);

  public:
    static constexpr result_type min() { return min_; }
    static constexpr result_type max() { return max_; }

    // Parameter and status classes
    class status_type {
      result_type r[int_math::ceil2(B)]{};
      unsigned int index{0};
      static constexpr unsigned int size() { return int_math::ceil2(B); }

    public:
      friend class lagfib2xor;

      // Equality comparable concept
      friend bool operator==(const status_type &a, const status_type &b) {
        if (a.index != b.index)
          return false;
        for (unsigned int i{0}; i < a.size(); ++i)
          if (a.r[i] != b.r[i])
            return false;
        return true;
      }

      friend bool operator!=(const status_type &a, const status_type &b) { return not(a == b); }

      // Streamable concept
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> &operator<<(
          std::basic_ostream<char_t, traits_t> &out, const status_type &S) {
        std::ios_base::fmtflags flags(out.flags());
        out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        out << '(' << S.index;
        for (unsigned int i{0}; i < S.size(); ++i)
          out << ' ' << S.r[i];
        out << ')';
        out.flags(flags);
        return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> &operator>>(
          std::basic_istream<char_t, traits_t> &in, status_type &S) {
        status_type S_new;
        std::ios_base::fmtflags flags(in.flags());
        in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
        in >> utility::delim('(') >> S_new.index;
        for (unsigned int i{0}; i < S.size(); ++i)
          in >> utility::delim(' ') >> S_new.r[i];
        in >> utility::delim(')');
        if (in)
          S = S_new;
        in.flags(flags);
        return in;
      }
    };

    // Random number engine concept
    lagfib2xor() { seed(); }

    explicit lagfib2xor(unsigned long s) { seed(s); }

    template<typename gen>
    explicit lagfib2xor(gen &g) {
      seed(g);
    }

    void seed() { seed(0); }

    void seed(unsigned long s) {
      minstd R(s);
      seed(R);
    }

    template<typename gen>
    void seed(gen &g) {
      for (unsigned int i{0}; i < B; ++i) {
        result_type r{0};
        for (int j{0}; j < std::numeric_limits<result_type>::digits; ++j) {
          r <<= 1;
          if (g() - gen::min() > gen::max() / 2)
            ++r;
        }
        S.r[i] = r;
      }
      S.index = B - 1;
    }

    void discard(unsigned long long n) {
      const unsigned int matrix_size = B;
      using matrix_type = matrix<GF2, matrix_size>;
      using vector_type = vector<result_type, matrix_size>;
      using size_type = typename matrix_type::size_type;
      const unsigned long long n_pivot{int_math::log2_ceil(n) * B * B * B};
      constexpr auto mask = int_math::mask(B);
      if (n > n_pivot) {
        const unsigned long long n_partial{n - matrix_size};
        matrix_type M;
        for (size_type i{0}; i < matrix_size - 1; ++i)
          M(i, i + 1) = GF2(true);
        M(matrix_size - 1, matrix_size - B) = GF2(true);
        M(matrix_size - 1, matrix_size - A) = GF2(true);
        M = power(M, n_partial);
        vector_type V;
        for (size_type i{0}; i < matrix_size; ++i)
          V(matrix_size - 1 - i) = S.r[(S.index - i) & mask];
        vector_type W;
        for (size_type i{0}; i < matrix_size; ++i) {
          result_type sum{};
          for (size_type k{0}; k < matrix_size; ++k)
            sum ^= M(i, k) * V(k);
          W(i) = sum;
        }
        S.index += n_partial;
        S.index &= mask;
        for (size_type i{0}; i < matrix_size; ++i)
          S.r[(S.index - i) & mask] = W(matrix_size - 1 - i);
        n -= n_partial;
      }
      for (unsigned long long i{0}; i < n; ++i)
        step();
    }

    // Equality comparable concept
    friend bool operator==(const lagfib2xor &R1, const lagfib2xor &R2) { return R1.S == R2.S; }

    friend bool operator!=(const lagfib2xor &R1, const lagfib2xor &R2) { return not(R1 == R2); }

    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const lagfib2xor &R) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '[' << lagfib2xor::name() << ' ' << R.S << ']';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, lagfib2xor &R) {
      typename lagfib2xor::status_type S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::ignore_spaces();
      in >> utility::delim('[') >> utility::delim(lagfib2xor::name()) >> utility::delim(' ') >>
          S_new >> utility::delim(']');
      if (in)
        R.S = S_new;
      in.flags(flags);
      return in;
    }

    // Other useful methods
  private:
    static std::string init_name() {
      std::stringstream name_str;
      name_str << "lagfib2xor_" << std::numeric_limits<result_type>::digits << '_' << A << '_'
               << B;
      return name_str.str();
    }

  public:
    static const char *name() {
      static const std::string name_str(init_name());
      return name_str.c_str();
    }
    long operator()(long x) {
      return static_cast<long>(utility::uniformco<double, lagfib2xor>(*this) * x);
    }

  private:
    status_type S;

    void step() {
      ++S.index;
      S.index &= int_math::mask(B);
      S.r[S.index] =
          S.r[(S.index - A) & int_math::mask(B)] ^ S.r[(S.index - B) & int_math::mask(B)];
    }
  };

  typedef lagfib2xor<unsigned long, 103, 250> r250_ul;
  typedef lagfib2xor<unsigned long long, 103, 250> r250_ull;
  typedef lagfib2xor<unsigned long, 168, 521> lagfib2xor_521_ul;
  typedef lagfib2xor<unsigned long long, 168, 521> lagfib2xor_521_ull;
  typedef lagfib2xor<unsigned long, 273, 607> lagfib2xor_607_ul;
  typedef lagfib2xor<unsigned long long, 273, 607> lagfib2xor_607_ull;
  typedef lagfib2xor<unsigned long, 418, 1279> lagfib2xor_1279_ul;
  typedef lagfib2xor<unsigned long long, 418, 1279> lagfib2xor_1279_ull;
  typedef lagfib2xor<unsigned long, 1029, 2281> lagfib2xor_2281_ul;
  typedef lagfib2xor<unsigned long long, 1029, 2281> lagfib2xor_2281_ull;
  typedef lagfib2xor<unsigned long, 576, 3217> lagfib2xor_3217_ul;
  typedef lagfib2xor<unsigned long long, 576, 3217> lagfib2xor_3217_ull;
  typedef lagfib2xor<unsigned long, 2098, 4423> lagfib2xor_4423_ul;
  typedef lagfib2xor<unsigned long long, 2098, 4423> lagfib2xor_4423_ull;
  typedef lagfib2xor<unsigned long, 4187, 9689> lagfib2xor_9689_ul;
  typedef lagfib2xor<unsigned long long, 4187, 9689> lagfib2xor_9689_ull;
  typedef lagfib2xor<unsigned long, 9842, 19937> lagfib2xor_19937_ul;
  typedef lagfib2xor<unsigned long long, 9842, 19937> lagfib2xor_19937_ull;

  typedef lagfib2xor<uint32_t, 103, 250> r250_32;
  typedef lagfib2xor<uint64_t, 103, 250> r250_64;
  typedef lagfib2xor<uint32_t, 168, 521> lagfib2xor_521_32;
  typedef lagfib2xor<uint64_t, 168, 521> lagfib2xor_521_64;
  typedef lagfib2xor<uint32_t, 273, 607> lagfib2xor_607_32;
  typedef lagfib2xor<uint64_t, 273, 607> lagfib2xor_607_64;
  typedef lagfib2xor<uint32_t, 418, 1279> lagfib2xor_1279_32;
  typedef lagfib2xor<uint64_t, 418, 1279> lagfib2xor_1279_64;
  typedef lagfib2xor<uint32_t, 1029, 2281> lagfib2xor_2281_32;
  typedef lagfib2xor<uint64_t, 1029, 2281> lagfib2xor_2281_64;
  typedef lagfib2xor<uint32_t, 576, 3217> lagfib2xor_3217_32;
  typedef lagfib2xor<uint64_t, 576, 3217> lagfib2xor_3217_64;
  typedef lagfib2xor<uint32_t, 2098, 4423> lagfib2xor_4423_32;
  typedef lagfib2xor<uint64_t, 2098, 4423> lagfib2xor_4423_64;
  typedef lagfib2xor<uint32_t, 4187, 9689> lagfib2xor_9689_32;
  typedef lagfib2xor<uint64_t, 4187, 9689> lagfib2xor_9689_64;
  typedef lagfib2xor<uint32_t, 9842, 19937> lagfib2xor_19937_32;
  typedef lagfib2xor<uint64_t, 9842, 19937> lagfib2xor_19937_64;

}  // namespace trng

#endif
