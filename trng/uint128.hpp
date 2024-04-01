// Copyright (c) 2000-2022, Heiko Bauke
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

#if !(defined TRNG_UINT128_HPP)

#define TRNG_UINT128_HPP

#include <cstdint>
#include <ciso646>
#include <ostream>
#include <istream>
#include <trng/int_types.hpp>


namespace trng {

  class uint128 {
#if defined __SIZEOF_INT128__
    __extension__ using uint128_t = unsigned __int128;
    uint128_t value{0};

    static uint128 create(uint128_t new_value) {
      uint128 result;
      result.value = new_value;
      return result;
    }
#else
    std::uint64_t m_lo{0};
    std::uint64_t m_hi{0};

    static uint128 create(std::uint64_t hi, std::uint64_t lo) { return uint128{hi, lo}; }
#endif

  public:
    uint128() = default;

#if defined __SIZEOF_INT128__
    explicit uint128(std::uint64_t lo) { value = lo; }
#else
    explicit uint128(std::uint64_t lo) : m_lo{lo} {}
#endif

#if defined __SIZEOF_INT128__
    explicit uint128(std::uint64_t hi, std::uint64_t lo) {
      value = hi;
      value <<= 64;
      value |= lo;
    }
#else
    explicit uint128(std::uint64_t hi, std::uint64_t lo) : m_lo{lo}, m_hi{hi} {}
#endif

    std::uint64_t lo() const {
#if defined __SIZEOF_INT128__
      return static_cast<std::uint64_t>(value & 0xffffffffffffffff);
#else
      return m_lo;
#endif
    }

    std::uint64_t hi() const {
#if defined __SIZEOF_INT128__
      return static_cast<std::uint64_t>(value >> 64);
#else
      return m_lo;
#endif
    }

    explicit operator std::uint64_t() const {
#if defined __SIZEOF_INT128__
      return static_cast<std::uint64_t>(value);
#else
      return m_lo;
#endif
    }

    explicit operator float() const {
#if defined __SIZEOF_INT128__
      return static_cast<float>(value);
#else
      return m_lo + 18446744073709551616.f * m_hi;
#endif
    }

    explicit operator double() const {
#if defined __SIZEOF_INT128__
      return static_cast<double>(value);
#else
      return m_lo + 18446744073709551616. * m_hi;
#endif
    }

    explicit operator long double() const {
#if defined __SIZEOF_INT128__
      return static_cast<long double>(value);
#else
      return m_lo + 18446744073709551616.l * m_hi;
#endif
    }

    friend inline bool operator==(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value == b.value;
#else
      return a.m_lo == b.m_lo and a.m_hi == b.m_hi;
#endif
    }

    friend inline bool operator!=(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value != b.value;
#else
      return a.m_lo != b.m_lo or a.m_hi != b.m_hi;
#endif
    }

    friend inline bool operator<(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value < b.value;
#else
      return a.m_hi < b.m_hi or (a.m_hi == b.m_hi and a.m_lo < b.m_lo);
#endif
    }

    friend inline bool operator<=(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value <= b.value;
#else
      bool y1 = a.m_hi < b.m_hi;
      bool y2 = a.m_hi == b.m_hi and a.m_lo <= b.m_lo;
      return y1 or y2;
      // return a.hi < b.hi or (a.hi == b.hi and a.lo <= b.lo);
#endif
    }

    friend inline bool operator>(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value > b.value;
#else
      return a.m_hi > b.m_hi or (a.m_hi == b.m_hi and a.m_lo > b.m_lo);
#endif
    }

    friend inline bool operator>=(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return a.value >= b.value;
#else
      return a.m_hi > b.m_hi or (a.m_hi == b.m_hi and a.m_lo >= b.m_lo);
#endif
    }

    uint128 operator+() const { return *this; }

    uint128 operator-() const {
#if defined __SIZEOF_INT128__
      return create(-value);
#else
      return ++create(~m_hi, ~m_lo);
#endif
    }

    uint128& operator+=(uint128 other) {
#if defined __SIZEOF_INT128__
      value += other.value;
      return *this;
#else
      m_lo += other.m_lo;
      m_hi += (m_lo < other.m_lo);
      m_hi += other.m_hi;
      return *this;
#endif
    }

    uint128& operator-=(uint128 other) {
#if defined __SIZEOF_INT128__
      value -= other.value;
      return *this;
#else
      const auto backup{m_lo};
      m_lo -= other.m_lo;
      m_hi -= (m_lo > backup);
      m_hi -= other.m_hi;
      return *this;
#endif
    }

    uint128& operator*=(uint128 other) {
#if defined __SIZEOF_INT128__
      value *= other.value;
      return *this;
#else
      const auto lo0{m_lo & 0xffffffff};
      const auto lo1{m_lo >> 32};
      const auto hi0{m_hi & 0xffffffff};
      const auto hi1{m_hi >> 32};

      const auto other_lo0{other.m_lo & 0xffffffff};
      const auto other_lo1{other.m_lo >> 32};
      const auto other_hi0{other.m_hi & 0xffffffff};
      const auto other_hi1{other.m_hi >> 32};

      const uint128 m0{0, lo0 * other_lo0};
      *this = m0;
      const uint128 m1{(lo1 * other_lo0) >> 32, (lo1 * other_lo0) << 32};
      *this += m1;
      const uint128 m2{(lo0 * other_lo1) >> 32, (lo0 * other_lo1) << 32};
      *this += m2;
      const std::uint64_t m3{hi0 * other_lo0};
      this->m_hi += m3;
      const std::uint64_t m4{lo0 * other_hi0};
      this->m_hi += m4;
      const std::uint64_t m5{lo1 * other_lo1};
      this->m_hi += m5;
      const std::uint64_t m6{(hi1 * other_lo0) << 32};
      this->m_hi += m6;
      const std::uint64_t m7{(lo0 * other_hi1) << 32};
      this->m_hi += m7;
      const std::uint64_t m8{(hi0 * other_lo1) << 32};
      this->m_hi += m8;
      const std::uint64_t m9{(lo1 * other_hi0) << 32};
      this->m_hi += m9;
      return *this;
#endif
    }

    uint128& operator/=(uint128 other) {
#if defined __SIZEOF_INT128__
      value /= other.value;
      return *this;
#else
      uint128 q;
      uint128 r;
      for (int i{127}; i >= 64; --i) {
        r <<= 1;
        r.m_lo |= (this->m_hi & (1ull << (i - 64))) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.m_hi |= 1ull << (i - 64);
        }
      }
      for (int i{63}; i >= 0; --i) {
        r <<= 1;
        r.m_lo |= (this->m_lo & (1ull << i)) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.m_lo |= 1ull << i;
        }
      }
      this->m_lo = q.m_lo;
      this->m_hi = q.m_hi;
      return *this;
#endif
    }

    uint128& operator%=(uint128 other) {
#if defined __SIZEOF_INT128__
      value %= other.value;
      return *this;
#else
      uint128 q;
      uint128 r;
      for (int i{127}; i >= 64; --i) {
        r <<= 1;
        r.m_lo |= (this->m_hi & (1ull << (i - 64))) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.m_hi |= 1ull << (i - 64);
        }
      }
      for (int i{63}; i >= 0; --i) {
        r <<= 1;
        r.m_lo |= (this->m_lo & (1ull << i)) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.m_lo |= 1ull << i;
        }
      }
      this->m_lo = r.m_lo;
      this->m_hi = r.m_hi;
      return *this;
#endif
    }

    friend inline uint128 operator+(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value + b.value);
#else
      return a += b;
#endif
    }

    friend inline uint128 operator-(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value - b.value);
#else
      return a -= b;
#endif
    }

    friend inline uint128 operator*(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value * b.value);
#else
      return a *= b;
#endif
    }

    friend inline uint128 operator/(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value / b.value);
#else
      return a /= b;
#endif
    }

    friend inline uint128 operator%(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value % b.value);
#else
      return a %= b;
#endif
    }

    uint128& operator++() {
#if defined __SIZEOF_INT128__
      ++value;
      return *this;
#else
      ++m_lo;
      if (m_lo == 0)
        ++m_hi;
      return *this;
#endif
    }

    uint128 operator++(int) {
      auto copy{*this};
      ++(*this);
      return copy;
    }

    uint128& operator--() {
#if defined __SIZEOF_INT128__
      --value;
      return *this;
#else
      if (m_lo == 0)
        --m_hi;
      --m_lo;
      return *this;
#endif
    }

    uint128 operator--(int) {
      auto copy{*this};
      --(*this);
      return copy;
    }

    uint128& operator|=(uint128 other) {
#if defined __SIZEOF_INT128__
      value |= other.value;
      return *this;
#else
      m_lo |= other.m_lo;
      m_hi |= other.m_hi;
      return *this;
#endif
    }

    uint128& operator&=(uint128 other) {
#if defined __SIZEOF_INT128__
      value &= other.value;
      return *this;
#else
      m_lo &= other.m_lo;
      m_hi &= other.m_hi;
      return *this;
#endif
    }

    uint128& operator^=(uint128 other) {
#if defined __SIZEOF_INT128__
      value ^= other.value;
      return *this;
#else
      m_lo ^= other.m_lo;
      m_hi ^= other.m_hi;
      return *this;
#endif
    }

    friend inline uint128 operator|(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value | b.value);
#else
      return create(a.m_hi | b.m_hi, a.m_lo | b.m_lo);
#endif
    }

    friend inline uint128 operator&(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value & b.value);
#else
      return create(a.m_hi & b.m_hi, a.m_lo & b.m_lo);
#endif
    }

    friend inline uint128 operator^(uint128 a, uint128 b) {
#if defined __SIZEOF_INT128__
      return create(a.value ^ b.value);
#else
      return create(a.m_hi ^ b.m_hi, a.m_lo ^ b.m_lo);
#endif
    }

    uint128 operator~() const {
#if defined __SIZEOF_INT128__
      return create(~value);
#else
      return create(~m_hi, ~m_lo);
#endif
    }

    uint128& operator>>=(int s) {
#if defined __SIZEOF_INT128__
      if (s < 0)
        return *this <<= (-s);
      value >>= s;
      return *this;
#else
      const uint128 result(*this >> s);
      *this = result;
      return *this;
#endif
    }

    uint128& operator<<=(int s) {
#if defined __SIZEOF_INT128__
      if (s < 0)
        return *this >>= (-s);
      value <<= s;
      return *this;
#else
      const uint128 result(*this << s);
      *this = result;
      return *this;
#endif
    }

    friend inline uint128 operator>>(uint128 a, int s) {
      if (s < 0)
        return a << (-s);
#if defined __SIZEOF_INT128__
      return create(a.value >> s);
#else
      if (s == 0)
        return a;
      if (s < 64)
        return uint128::create(a.m_hi >> s, (a.m_lo >> s) | (a.m_hi << (64 - s)));
      return uint128::create(0, a.m_hi >> (s - 64));
#endif
    }

    friend inline uint128 operator<<(uint128 a, int s) {
      if (s < 0)
        return a >> (-s);
#if defined __SIZEOF_INT128__
      return create(a.value << s);
#else
      if (s == 0)
        return a;
      if (s < 64)
        return uint128::create((a.m_hi << s) | (a.m_lo >> (64 - s)), a.m_lo << s);
      return uint128::create(a.m_lo << (s - 64), 0);
#endif
    }
  };


  template<typename T, typename TRAIT>
  std::basic_ostream<T, TRAIT>& operator<<(std::basic_ostream<T, TRAIT>& out, uint128 value) {
    const T zero_c{static_cast<T>('0')};
    const uint128 zero{0};
    const uint128 ten{10};
    const uint128 seven{7};
    const uint128 fifteen{15};

    const auto to_hex = [](std::uint64_t x, bool upper_case) -> T {
      if (x < 10)
        return static_cast<T>('0' + x);
      else {
        if (upper_case)
          return static_cast<T>('A' + x - 10);
        else
          return static_cast<T>('a' + x - 10);
      }
    };

    const auto stream_flags{out.flags()};
    const bool upper_case{(stream_flags & std::ios_base::uppercase) ==
                          std::ios_base::uppercase};
    const bool showbase{(stream_flags & std::ios_base::showbase) == std::ios_base::showbase};

    std::basic_string<T> digits;
    std::basic_string<T> padding;
    std::basic_string<T> base;

    if (value == zero) {
      digits.push_back(static_cast<T>('0'));
    } else {
      switch (stream_flags & std::ios_base::basefield) {
        case std::ios_base::oct:
          while (value > zero) {
            digits.insert(0, 1, to_hex(static_cast<std::uint64_t>(value & seven), upper_case));
            value >>= 3;
          }
          break;
        case std::ios_base::dec:
        default:
          while (value > zero) {
            digits.insert(0, 1,
                          static_cast<T>(static_cast<std::uint64_t>(value % ten)) + zero_c);
            value /= ten;
          }
          break;
        case std::ios_base::hex:
          while (value > zero) {
            digits.insert(0, 1,
                          to_hex(static_cast<std::uint64_t>(value & fifteen), upper_case));
            value >>= 4;
          }
          break;
      }
    }

    if (showbase) {
      switch (stream_flags & std::ios_base::basefield) {
        case std::ios_base::oct:
          base.push_back(static_cast<T>('0'));
          break;
        case std::ios_base::dec:
        default:
          break;
        case std::ios_base::hex:
          base.push_back(static_cast<T>('0'));
          base.push_back(static_cast<T>(upper_case ? 'X' : 'x'));
          break;
      }
    }

    if (static_cast<std::streamsize>(base.size() + digits.size()) < out.width())
      padding.insert(0, out.width() - base.size() - digits.size(), out.fill());

    switch (stream_flags & std::ios_base::adjustfield) {
      case std::ios_base::left:
      default:
        out.write(base.c_str(), base.size());
        out.write(digits.c_str(), digits.size());
        out.write(padding.c_str(), padding.size());
        break;
      case std::ios_base::internal:
        out.write(base.c_str(), base.size());
        out.write(padding.c_str(), padding.size());
        out.write(digits.c_str(), digits.size());
        break;
      case std::ios_base::right:
        out.write(padding.c_str(), padding.size());
        out.write(base.c_str(), base.size());
        out.write(digits.c_str(), digits.size());
        break;
    }
    return out;
  }


  template<typename T, typename TRAIT>
  std::basic_istream<T, TRAIT>& operator>>(std::basic_istream<T, TRAIT>& in, uint128& value) {
    const auto stream_flags{in.flags()};

    const auto is_dec_digit = [](typename TRAIT::int_type c) -> bool {
      return T('0') <= c and c <= T('9');
    };
    const auto dec_digit_to_uint128 = [](typename TRAIT::int_type c) -> uint128 {
      return static_cast<uint128>(
          static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
    };

    const auto is_oct_digit = [](typename TRAIT::int_type c) -> bool {
      return T('0') <= c and c <= T('7');
    };
    const auto oct_digit_to_uint128 = [](typename TRAIT::int_type c) -> uint128 {
      return static_cast<uint128>(
          static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
    };

    const auto is_hex_digit = [](typename TRAIT::int_type c) -> bool {
      return (T('0') <= c and c <= T('9')) or (T('a') <= c and c <= T('f')) or
             (T('A') <= c and c <= T('F'));
    };
    const auto is_hex_basechar = [](typename TRAIT::int_type c) -> bool {
      return c == T('x') or c == T('X');
    };
    const auto hex_digit_to_uint128 = [](typename TRAIT::int_type c) -> uint128 {
      if (T('0') <= c and c <= T('9'))
        return static_cast<uint128>(
            static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
      if (T('a') <= c and c <= T('f'))
        return static_cast<uint128>(
            static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('a') + 10));
      return static_cast<uint128>(
          static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('A') + 10));
    };

    in >> std::ws;

    switch (stream_flags & std::ios_base::basefield) {
      case std::ios_base::oct: {
        const uint128 mult_overflow{0x2000000000000000, 0x000000000000000};  // 2^128 / 8
        if (not is_oct_digit(in.peek())) {
          in.setstate(std::ios_base::failbit);
          return in;
        }
        uint128 result;
        while (is_oct_digit(in.peek())) {
          if (result >= mult_overflow) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          result *= static_cast<uint128>(static_cast<std::uint64_t>(8));
          result += oct_digit_to_uint128(in.get());
        }
        value = result;
      } break;
      case std::ios_base::dec:
      default: {
        const uint128 mult_overflow{0x1999999999999999, 0x999999999999999a};  // ceil(2^128 / 10)
        const uint128 max{0xffffffffffffffff, 0xffffffffffffffff};  // 2^128 - 1
        if (not is_dec_digit(in.peek())) {
          in.setstate(std::ios_base::failbit);
          return in;
        }
        uint128 result;
        while (is_dec_digit(in.peek())) {
          if (result >= mult_overflow) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          result *= static_cast<uint128>(static_cast<std::uint64_t>(10));
          const auto next_digit{dec_digit_to_uint128(in.get())};
          if (next_digit > max - result) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          result += next_digit;
        }
        value = result;
      } break;
      case std::ios_base::hex: {
        const uint128 mult_overflow{0x1000000000000000, 0x000000000000000};  // 2^128 / 16
        if (not is_hex_digit(in.peek())) {
          in.setstate(std::ios_base::failbit);
          return in;
        }
        uint128 result;
        int pos{0};
        while (is_hex_digit(in.peek()) or (pos == 1 and is_hex_basechar(in.peek()))) {
          if (pos == 1 and is_hex_basechar(in.peek())) {
            in.get();
          } else {
            if (result >= mult_overflow) {
              in.setstate(std::ios_base::failbit);
              return in;
            }
            result *= static_cast<uint128>(static_cast<std::uint64_t>(16));
            result += hex_digit_to_uint128(in.get());
          }
          ++pos;
        }
        value = result;
      } break;
    }

    return in;
  }

}  // namespace trng

#endif
