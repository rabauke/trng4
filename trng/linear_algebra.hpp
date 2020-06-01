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

#if !(defined TRNG_LINEAR_ALGEBRA_HPP)

#define TRNG_LINEAR_ALGEBRA_HPP

#include <vector>
#include <cstddef>
#include <ciso646>
#include <iostream>
#include <type_traits>
#include <trng/utility.hpp>
#include <trng/int_math.hpp>

namespace trng {

  template<typename T, std::size_t n>
  class vector {
    std::vector<T> data;

  public:
    using size_type = typename std::vector<T>::size_type;
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    vector() : data(n) {}

    template<typename F>
    explicit vector(F f) {
      static_assert(trng::utility::is_same<T, decltype(f(0))>::value,
                    "wrong return type of functor");
      data.reserve(n);
      for (size_type i{0}; i < n; ++i)
        data.push_back(f(i));
    }

    template<typename... Ts>
    explicit vector(Ts... t) : data{t...} {
      static_assert(sizeof...(Ts) == n, "wrong number of arguments");
      static_assert(trng::utility::is_same<T, Ts...>::value,
                    "wrong type in constructor argument");
    }

    reference operator()(size_type i) { return data[i]; }
    const_reference operator()(size_type i) const { return data[i]; }
    constexpr size_type size() const { return n; }
    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }
    iterator end() { return data.end(); }
    const_iterator end() const { return data.end(); }
    bool operator==(const vector &other) const { return data == other.data; }
    bool operator!=(const vector &other) const { return data != other.data; }
  };


  template<typename T, std::size_t n>
  class matrix {
    std::vector<T> data;

  public:
    using size_type = typename std::vector<T>::size_type;
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;

    matrix() : data(n * n) {}

    template<typename F>
    explicit matrix(F f) {
      static_assert(trng::utility::is_same<T, decltype(f(0, 0))>::value,
                    "wrong return type of functor");
      data.reserve(n * n);
      for (size_type i{0}; i < n; ++i)
        for (size_type j{0}; j < n; ++j)
          data.push_back(f(i, j));
    }

    template<typename... Ts>
    explicit matrix(Ts... t) : data{t...} {
      static_assert(sizeof...(Ts) == n * n, "wrong number of arguments");
      static_assert(trng::utility::is_same<T, Ts...>::value,
                    "wrong type in constructor argument");
    }

    reference operator()(size_type i, size_type j) { return data[j + i * n]; }
    const_reference operator()(size_type i, size_type j) const { return data[j + i * n]; }
    constexpr size_type size() const { return n; }
    bool operator==(const matrix &other) const { return data == other.data; }
    bool operator!=(const matrix &other) const { return data != other.data; }
  };


  template<typename T, std::size_t n>
  vector<T, n> operator+(const vector<T, n> &a, const vector<T, n> &b) {
    using size_type = typename vector<T, n>::size_type;
    auto f = [&](size_type i) -> T { return a(i) + b(i); };
    return vector<T, n>(f);
  }


  template<typename T, std::size_t n>
  vector<T, n> operator*(T a, const vector<T, n> &b) {
    using size_type = typename vector<T, n>::size_type;
    auto f = [&](size_type i) -> T { return a * b(i); };
    return vector<T, n>(f);
  }


  template<typename T, std::size_t n>
  vector<T, n> operator*(const vector<T, n> &a, T b) {
    using size_type = typename vector<T, n>::size_type;
    auto f = [&](size_type i) -> T { return a(i) * b; };
    return vector<T, n>(f);
  }


  template<typename T, std::size_t n>
  matrix<T, n> operator+(const matrix<T, n> &a, const matrix<T, n> &b) {
    using size_type = typename matrix<T, n>::size_type;
    auto f = [&](size_type i, size_type j) -> T { return a(i, j) + b(i, j); };
    return matrix<T, n>(f);
  }


  template<typename T, std::size_t n>
  matrix<T, n> operator*(T a, const matrix<T, n> &b) {
    using size_type = typename matrix<T, n>::size_type;
    auto f = [&](size_type i, size_type j) -> T { return a * b(i, j); };
    return matrix<T, n>(f);
  }


  template<typename T, std::size_t n>
  matrix<T, n> operator*(const matrix<T, n> &a, T b) {
    using size_type = typename matrix<T, n>::size_type;
    auto f = [&](size_type i, size_type j) -> T { return a(i, j) * b; };
    return matrix<T, n>(f);
  }


  template<typename T, std::size_t n>
  vector<T, n> operator*(const matrix<T, n> &a, const vector<T, n> &b) {
    using size_type = typename matrix<T, n>::size_type;
    auto f = [&](size_type i) -> T {
      T sum{};
      for (size_type k{0}; k < n; ++k)
        sum = sum + a(i, k) * b(k);
      return sum;
    };
    return vector<T, n>(f);
  }


  template<typename T, std::size_t n>
  matrix<T, n> operator*(const matrix<T, n> &a, const matrix<T, n> &b) {
    using size_type = typename matrix<T, n>::size_type;
    matrix<T, n> res;
    const size_type step{32};
    for (size_type j0{0}; j0 < n; j0 += step) {
      for (size_type i{0}; i < n; ++i)
        for (size_type j{j0}, j_end{utility::min(n, j0 + step)}; j < j_end; ++j)
          res(i, j) = T{0};
      for (size_type k0{0}; k0 < n; k0 += step) {
        for (size_type i{0}; i < n; ++i)
          for (size_type j{j0}, j_end{utility::min(n, j0 + step)}; j < j_end; ++j) {
            T sum{0};
            for (size_type k{k0}, k_end{utility::min(n, k0 + step)}; k < k_end; ++k)
              sum = sum + a(i, k) * b(k, j);
            res(i, j) = res(i, j) + sum;
          }
      }
    }
    return res;
  }


  template<typename T, std::size_t n>
  matrix<T, n> power(const matrix<T, n> &a, unsigned long long m) {
    using matrix_type = matrix<T, n>;
    using size_type = typename matrix_type::size_type;

    auto unit = [&](size_type i, size_type j) -> T { return i == j ? T{1} : T{0}; };
    matrix_type res(unit);
    matrix_type powers(a);
    while (m > 0) {
      if ((m & 1ull) == 1ull)
        res = res * powers;
      m >>= 1u;
      if (m == 0)
        break;
      powers = powers * powers;
    }
    return res;
  }


  class GF2 {
  public:
    using value_type = std::uint8_t;

  private:
    std::uint8_t value{0};

    explicit GF2(value_type v) : value{v} {}

  public:
    explicit GF2(bool v = false) : value(v ? 1 : 0) {}
    explicit GF2(int v) : value(v != 0 ? 1 : 0) {}

    explicit operator bool() { return value; }

    friend bool operator==(const GF2 a, const GF2 b) { return a.value == b.value; }
    friend bool operator!=(const GF2 a, const GF2 b) { return a.value != b.value; }

    GF2 &operator+=(const GF2 other) {
      value ^= other.value;
      return *this;
    }

    GF2 &operator*=(const GF2 other) {
      value &= other.value;
      return *this;
    }

    friend GF2 operator+(const GF2 a, const GF2 b) {
      return GF2{static_cast<value_type>(a.value ^ b.value)};
    }

    friend GF2 operator*(const GF2 a, const GF2 b) {
      return GF2{static_cast<value_type>(a.value & b.value)};
    }

    template<typename T>
    friend typename std::enable_if<std::is_integral<T>::value, T>::type operator*(const GF2 a,
                                                                                  const T b) {
      return a.value ? b : T{};
    }

    template<typename T>
    friend typename std::enable_if<std::is_integral<T>::value, T>::type operator*(const T a,
                                                                                  const GF2 b) {
      return b.value ? a : T{};
    }
  };


  template<typename CharT, typename Traits>
  std::basic_ostream<CharT, Traits> &operator<<(std::basic_ostream<CharT, Traits> &os,
                                                GF2 value) {
    return os << static_cast<bool>(value);
  }


}  // namespace trng

#endif
