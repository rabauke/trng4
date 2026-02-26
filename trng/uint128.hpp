// Copyright (c) 2000-2026, Heiko Bauke
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
#if defined _MSC_VER && __cplusplus <= 201703
#include <ciso646>
#endif
#include <ostream>
#include <istream>
#include <trng/int_types.hpp>


namespace trng {

  namespace impl {

    template<typename T, typename TRAIT, typename UINT128>
    std::basic_ostream<T,
        TRAIT> &operator<<(std::basic_ostream<T,
        TRAIT> &out, UINT128 value) {
      const T zero_c{static_cast<T>('0')};
      const UINT128 zero{0};
      const UINT128 ten{10};
      const UINT128 seven{7};
      const UINT128 fifteen{15};

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
              digits.insert(0, 1,
                            to_hex(static_cast<std::uint64_t>(value & seven), upper_case));
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


    template<typename T, typename TRAIT, typename UINT128>
    std::basic_istream<T,
        TRAIT> &operator>>(std::basic_istream<T,
        TRAIT> &in, UINT128 &value) {
      const auto stream_flags{in.flags()};

      const auto is_dec_digit = [](typename TRAIT::int_type c) -> bool {
        return T('0') <= c and c <= T('9');
      };
      const auto dec_digit_to_uint128 = [](typename TRAIT::int_type c) -> UINT128 {
        return static_cast<UINT128>(
            static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
      };

      const auto is_oct_digit = [](typename TRAIT::int_type c) -> bool {
        return T('0') <= c and c <= T('7');
      };
      const auto oct_digit_to_uint128 = [](typename TRAIT::int_type c) -> UINT128 {
        return static_cast<UINT128>(
            static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
      };

      const auto is_hex_digit = [](typename TRAIT::int_type c) -> bool {
        return (T('0') <= c and c <= T('9')) or (T('a') <= c and c <= T('f')) or
               (T('A') <= c and c <= T('F'));
      };
      const auto is_hex_basechar = [](typename TRAIT::int_type c) -> bool {
        return c == T('x') or c == T('X');
      };
      const auto hex_digit_to_uint128 = [](typename TRAIT::int_type c) -> UINT128 {
        if (T('0') <= c and c <= T('9'))
          return static_cast<UINT128>(
              static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('0')));
        if (T('a') <= c and c <= T('f'))
          return static_cast<UINT128>(
              static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('a') + 10));
        return static_cast<UINT128>(
            static_cast<std::uint64_t>(c - static_cast<typename TRAIT::int_type>('A') + 10));
      };

      in >> std::ws;

      switch (stream_flags & std::ios_base::basefield) {
        case std::ios_base::oct: {
          const UINT128 mult_overflow{0x2000000000000000, 0x000000000000000};  // 2^128 / 8
          if (not is_oct_digit(in.peek())) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          UINT128 result;
          while (is_oct_digit(in.peek())) {
            if (result >= mult_overflow) {
              in.setstate(std::ios_base::failbit);
              return in;
            }
            result *= static_cast<UINT128>(static_cast<std::uint64_t>(8));
            result += oct_digit_to_uint128(in.get());
          }
          value = result;
        }
          break;
        case std::ios_base::dec:
        default: {
          const UINT128 mult_overflow{0x1999999999999999,
                                      0x999999999999999a};  // ceil(2^128 / 10)
          const UINT128 max{0xffffffffffffffff, 0xffffffffffffffff};  // 2^128 - 1
          if (not is_dec_digit(in.peek())) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          UINT128 result;
          while (is_dec_digit(in.peek())) {
            if (result >= mult_overflow) {
              in.setstate(std::ios_base::failbit);
              return in;
            }
            result *= static_cast<UINT128>(static_cast<std::uint64_t>(10));
            const auto next_digit{dec_digit_to_uint128(in.get())};
            if (next_digit > max - result) {
              in.setstate(std::ios_base::failbit);
              return in;
            }
            result += next_digit;
          }
          value = result;
        }
          break;
        case std::ios_base::hex: {
          const UINT128 mult_overflow{0x1000000000000000, 0x000000000000000};  // 2^128 / 16
          if (not is_hex_digit(in.peek())) {
            in.setstate(std::ios_base::failbit);
            return in;
          }
          UINT128 result;
          int pos{0};
          while (is_hex_digit(in.peek()) or (pos == 1 and is_hex_basechar(in.peek()))) {
            if (pos == 1 and is_hex_basechar(in.peek())) {
              in.get();
            } else {
              if (result >= mult_overflow) {
                in.setstate(std::ios_base::failbit);
                return in;
              }
              result *= static_cast<UINT128>(static_cast<std::uint64_t>(16));
              result += hex_digit_to_uint128(in.get());
            }
            ++pos;
          }
          value = result;
        }
          break;
      }

      return in;
    }

  }


  namespace optimized_impl {

#if defined __SIZEOF_INT128__

    class uint128 {
      __extension__ using uint128_t = unsigned __int128;
      uint128_t value{0};

      static uint128 create(uint128_t new_value) {
        uint128 result;
        result.value = new_value;
        return result;
      }

    public:
      uint128() = default;

      explicit uint128(std::uint64_t lo) { value = lo; }

      explicit uint128(std::uint64_t hi, std::uint64_t lo) {
        value = hi;
        value <<= 64;
        value |= lo;
      }

      std::uint64_t lo() const {
        return static_cast<std::uint64_t>(value & 0xffffffffffffffff);
      }

      std::uint64_t hi() const {
        return static_cast<std::uint64_t>(value >> 64);
      }

      explicit operator std::uint64_t() const {
        return static_cast<std::uint64_t>(value);
      }

      explicit operator float() const {
        return static_cast<float>(value);
      }

      explicit operator double() const {
        return static_cast<double>(value);
      }

      explicit operator long double() const {
        return static_cast<long double>(value);
      }

      friend inline bool operator==(uint128 a, uint128 b) {
        return a.value == b.value;
      }

      friend inline bool operator!=(uint128 a, uint128 b) {
        return a.value != b.value;
      }

      friend inline bool operator<(uint128 a, uint128 b) {
        return a.value < b.value;
      }

      friend inline bool operator<=(uint128 a, uint128 b) {
        return a.value <= b.value;
      }

      friend inline bool operator>(uint128 a, uint128 b) {
        return a.value > b.value;
      }

      friend inline bool operator>=(uint128 a, uint128 b) {
        return a.value >= b.value;
      }

      uint128 operator+() const { return *this; }

      uint128 operator-() const {
        return create(-value);
      }

      uint128 &operator+=(uint128 other) {
        value += other.value;
        return *this;
      }

      uint128 &operator-=(uint128 other) {
        value -= other.value;
        return *this;
      }

      uint128 &operator*=(uint128 other) {
        value *= other.value;
        return *this;
      }

      uint128 &operator/=(uint128 other) {
        value /= other.value;
        return *this;
      }

      uint128 &operator%=(uint128 other) {
        value %= other.value;
        return *this;
      }

      friend inline uint128 operator+(uint128 a, uint128 b) {
        return create(a.value + b.value);
      }

      friend inline uint128 operator-(uint128 a, uint128 b) {
        return create(a.value - b.value);
      }

      friend inline uint128 operator*(uint128 a, uint128 b) {
        return create(a.value * b.value);
      }

      friend inline uint128 operator/(uint128 a, uint128 b) {
        return create(a.value / b.value);
      }

      friend inline uint128 operator%(uint128 a, uint128 b) {
        return create(a.value % b.value);
      }

      uint128 &operator++() {
        ++value;
        return *this;
      }

      uint128 operator++(int) {
        auto copy{*this};
        ++(*this);
        return copy;
      }

      uint128 &operator--() {
        --value;
        return *this;
      }

      uint128 operator--(int) {
        auto copy{*this};
        --(*this);
        return copy;
      }

      uint128 &operator|=(uint128 other) {
        value |= other.value;
        return *this;
      }

      uint128 &operator&=(uint128 other) {
        value &= other.value;
        return *this;
      }

      uint128 &operator^=(uint128 other) {
        value ^= other.value;
        return *this;
      }

      friend inline uint128 operator|(uint128 a, uint128 b) {
        return create(a.value | b.value);
      }

      friend inline uint128 operator&(uint128 a, uint128 b) {
        return create(a.value & b.value);
      }

      friend inline uint128 operator^(uint128 a, uint128 b) {
        return create(a.value ^ b.value);
      }

      uint128 operator~() const {
        return create(~value);
      }

      uint128 &operator>>=(int s) {
        if (s < 0)
          return *this <<= (-s);
        value >>= s;
        return *this;
      }

      uint128 &operator<<=(int s) {
        if (s < 0)
          return *this >>= (-s);
        value <<= s;
        return *this;
      }

      friend inline uint128 operator>>(uint128 a, int s) {
        if (s < 0)
          return a << (-s);
        return create(a.value >> s);
      }

      friend inline uint128 operator<<(uint128 a, int s) {
        if (s < 0)
          return a >> (-s);
        return create(a.value << s);
      }
    };


    template<typename T, typename TRAIT>
    std::basic_ostream<T,
        TRAIT> &operator<<(std::basic_ostream<T,
        TRAIT> &out, const uint128 &value) {
      return trng::impl::operator<<(out, value);
    }


    template<typename T, typename TRAIT>
    std::basic_istream<T,
        TRAIT> &operator>>(std::basic_istream<T,
        TRAIT> &in, uint128 &value) {
      return trng::impl::operator>>(in, value);
    }

#endif

  }


  namespace portable_impl {

    class uint128 {
      std::uint64_t m_lo{0};
      std::uint64_t m_hi{0};

      static uint128 create(std::uint64_t hi, std::uint64_t lo) { return uint128{hi, lo}; }

    public:
      uint128() = default;

      explicit uint128(std::uint64_t lo) :
          m_lo{lo} { }

      explicit uint128(std::uint64_t hi, std::uint64_t lo) :
          m_lo{lo}, m_hi{hi} { }

      std::uint64_t lo() const {
        return m_lo;
      }

      std::uint64_t hi() const {
        return m_lo;
      }

      explicit operator std::uint64_t() const {
        return m_lo;
      }

      explicit operator float() const {
        return m_lo + 18446744073709551616.f * m_hi;
      }

      explicit operator double() const {
        return m_lo + 18446744073709551616. * m_hi;
      }

      explicit operator long double() const {
        return m_lo + 18446744073709551616.l * m_hi;
      }

      friend inline bool operator==(uint128 a, uint128 b) {
        return a.m_lo == b.m_lo and a.m_hi == b.m_hi;
      }

      friend inline bool operator!=(uint128 a, uint128 b) {
        return a.m_lo != b.m_lo or a.m_hi != b.m_hi;
      }

      friend inline bool operator<(uint128 a, uint128 b) {
        return a.m_hi < b.m_hi or (a.m_hi == b.m_hi and a.m_lo < b.m_lo);
      }

      friend inline bool operator<=(uint128 a, uint128 b) {
        return a.m_hi < b.m_hi or (a.m_hi == b.m_hi and a.m_lo <= b.m_lo);
      }

      friend inline bool operator>(uint128 a, uint128 b) {
        return a.m_hi > b.m_hi or (a.m_hi == b.m_hi and a.m_lo > b.m_lo);
      }

      friend inline bool operator>=(uint128 a, uint128 b) {
        return a.m_hi > b.m_hi or (a.m_hi == b.m_hi and a.m_lo >= b.m_lo);
      }

      uint128 operator+() const { return *this; }

      uint128 operator-() const {
        return ++create(~m_hi, ~m_lo);
      }

      uint128 &operator+=(uint128 other) {
        m_lo += other.m_lo;
        m_hi += (m_lo < other.m_lo);
        m_hi += other.m_hi;
        return *this;
      }

      uint128 &operator-=(uint128 other) {
        const auto backup{m_lo};
        m_lo -= other.m_lo;
        m_hi -= (m_lo > backup);
        m_hi -= other.m_hi;
        return *this;
      }

      uint128 &operator*=(uint128 other) {
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
      }

      uint128 &operator/=(uint128 other) {
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
      }

      uint128 &operator%=(uint128 other) {
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
      }

      friend inline uint128 operator+(uint128 a, uint128 b) {
        return a += b;
      }

      friend inline uint128 operator-(uint128 a, uint128 b) {
        return a -= b;
      }

      friend inline uint128 operator*(uint128 a, uint128 b) {
        return a *= b;
      }

      friend inline uint128 operator/(uint128 a, uint128 b) {
        return a /= b;
      }

      friend inline uint128 operator%(uint128 a, uint128 b) {
        return a %= b;
      }

      uint128 &operator++() {
        ++m_lo;
        if (m_lo == 0)
          ++m_hi;
        return *this;
      }

      uint128 operator++(int) {
        auto copy{*this};
        ++(*this);
        return copy;
      }

      uint128 &operator--() {
        if (m_lo == 0)
          --m_hi;
        --m_lo;
        return *this;
      }

      uint128 operator--(int) {
        auto copy{*this};
        --(*this);
        return copy;
      }

      uint128 &operator|=(uint128 other) {
        m_lo |= other.m_lo;
        m_hi |= other.m_hi;
        return *this;
      }

      uint128 &operator&=(uint128 other) {
        m_lo &= other.m_lo;
        m_hi &= other.m_hi;
        return *this;
      }

      uint128 &operator^=(uint128 other) {
        m_lo ^= other.m_lo;
        m_hi ^= other.m_hi;
        return *this;
      }

      friend inline uint128 operator|(uint128 a, uint128 b) {
        return create(a.m_hi | b.m_hi, a.m_lo | b.m_lo);
      }

      friend inline uint128 operator&(uint128 a, uint128 b) {
        return create(a.m_hi & b.m_hi, a.m_lo & b.m_lo);
      }

      friend inline uint128 operator^(uint128 a, uint128 b) {
        return create(a.m_hi ^ b.m_hi, a.m_lo ^ b.m_lo);
      }

      uint128 operator~() const {
        return create(~m_hi, ~m_lo);
      }

      uint128 &operator>>=(int s) {
        const uint128 result(*this >> s);
        *this = result;
        return *this;
      }

      uint128 &operator<<=(int s) {
        const uint128 result(*this << s);
        *this = result;
        return *this;
      }

      friend inline uint128 operator>>(uint128 a, int s) {
        if (s < 0)
          return a << (-s);
        if (s == 0)
          return a;
        if (s < 64)
          return uint128::create(a.m_hi >> s, (a.m_lo >> s) | (a.m_hi << (64 - s)));
        return uint128::create(0, a.m_hi >> (s - 64));
      }

      friend inline uint128 operator<<(uint128 a, int s) {
        if (s < 0)
          return a >> (-s);
        if (s == 0)
          return a;
        if (s < 64)
          return uint128::create((a.m_hi << s) | (a.m_lo >> (64 - s)), a.m_lo << s);
        return uint128::create(a.m_lo << (s - 64), 0);
      }
    };


    template<typename T, typename TRAIT>
    std::basic_ostream<T,
        TRAIT> &operator<<(std::basic_ostream<T,
        TRAIT> &out, const uint128 &value) {
      return trng::impl::operator<<(out, value);
    }


    template<typename T, typename TRAIT>
    std::basic_istream<T,
        TRAIT> &operator>>(std::basic_istream<T,
        TRAIT> &in, uint128 &value) {
      return trng::impl::operator>>(in, value);
    }

  }


#if defined __SIZEOF_INT128__
  using namespace optimized_impl;
#else
  using namespace portable_impl;
#endif

}  // namespace trng

#endif
