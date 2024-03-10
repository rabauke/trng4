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
    uint64_t lo{0};
    uint64_t hi{0};

    static uint128 create(uint64_t hi, uint64_t lo) { return uint128{hi, lo}; }
#endif

  public:
    uint128() = default;

#if defined __SIZEOF_INT128__
    explicit uint128(uint64_t lo) { value = lo; }
#else
    explicit uint128(uint64_t lo) : lo{lo} {}
#endif

#if defined __SIZEOF_INT128__
    explicit uint128(uint64_t hi, uint64_t lo) {
      value = hi;
      value <<= 64;
      value |= lo;
    }
#else
    explicit uint128(uint64_t hi, uint64_t lo) : lo{lo}, hi{hi} {}
#endif

#if defined __SIZEOF_INT128__
    bool operator==(uint128 other) const { return value == other.value; }
#else
    bool operator==(uint128 other) const { return lo == other.lo and hi == other.hi; }
#endif

#if defined __SIZEOF_INT128__
    bool operator!=(uint128 other) const { return value != other.value; }
#else
    bool operator!=(uint128 other) const { return lo != other.lo or hi != other.hi; }
#endif

#if defined __SIZEOF_INT128__
    bool operator<(uint128 other) const { return value < other.value; }
#else
    bool operator<(uint128 other) const {
      return hi < other.hi or (hi == other.hi and lo < other.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    bool operator<=(uint128 other) const { return value <= other.value; }
#else
    bool operator<=(uint128 other) const {
      return hi < other.hi or (hi == other.hi and lo <= other.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    bool operator>(uint128 other) const { return value > other.value; }
#else
    bool operator>(uint128 other) const {
      return hi > other.hi or (hi == other.hi and lo > other.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    bool operator>=(uint128 other) const { return value >= other.value; }
#else
    bool operator>=(uint128 other) const {
      return hi > other.hi or (hi == other.hi and lo >= other.lo);
    }
#endif


#if defined __SIZEOF_INT128__
    uint128& operator+=(uint128 other) {
      value += other.value;
      return *this;
    }
#else
    uint128& operator+=(uint128 other) {
      lo += other.lo;
      hi += (lo < other.lo);
      hi += other.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator-=(uint128 other) {
      value -= other.value;
      return *this;
    }
#else
    uint128& operator-=(uint128 other) {
      const auto backup{lo};
      lo -= other.lo;
      hi -= (lo > backup);
      hi -= other.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator*=(uint128 other) {
      value *= other.value;
      return *this;
    }
#else
    uint128& operator*=(uint128 other) {
      const auto lo0{lo & 0xffffffff};
      const auto lo1{lo >> 32};
      const auto hi0{hi & 0xffffffff};
      const auto hi1{hi >> 32};

      const auto other_lo0{other.lo & 0xffffffff};
      const auto other_lo1{other.lo >> 32};
      const auto other_hi0{other.hi & 0xffffffff};
      const auto other_hi1{other.hi >> 32};

      const uint128 m0{0, lo0 * other_lo0};
      *this = m0;
      const uint128 m1{(lo1 * other_lo0) >> 32, (lo1 * other_lo0) << 32};
      *this += m1;
      const uint128 m2{(lo0 * other_lo1) >> 32, (lo0 * other_lo1) << 32};
      *this += m2;
      const u_int64_t m3{hi0 * other_lo0};
      this->hi += m3;
      const u_int64_t m4{lo0 * other_hi0};
      this->hi += m4;
      const u_int64_t m5{lo1 * other_lo1};
      this->hi += m5;
      const u_int64_t m6{(hi1 * other_lo0) << 32};
      this->hi += m6;
      const u_int64_t m7{(lo0 * other_hi1) << 32};
      this->hi += m7;
      const u_int64_t m8{(hi0 * other_lo1) << 32};
      this->hi += m8;
      const u_int64_t m9{(lo1 * other_hi0) << 32};
      this->hi += m9;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator/=(uint128 other) {
      value /= other.value;
      return *this;
    }
#else
    uint128& operator/=(uint128 other) {
      uint128 q;
      uint128 r;
      for (int i{127}; i >= 64; --i) {
        r <<= 1;
        r.lo |= (this->hi & (1ull << (i - 64))) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.hi |= 1ull << (i - 64);
        }
      }
      for (int i{63}; i >= 0; --i) {
        r <<= 1;
        r.lo |= (this->lo & (1ull << i)) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.lo |= 1ull << i;
        }
      }
      this->lo = q.lo;
      this->hi = q.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator%=(uint128 other) {
      value %= other.value;
      return *this;
    }
#else
    uint128& operator%=(uint128 other) {
      uint128 q;
      uint128 r;
      for (int i{127}; i >= 64; --i) {
        r <<= 1;
        r.lo |= (this->hi & (1ull << (i - 64))) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.hi |= 1ull << (i - 64);
        }
      }
      for (int i{63}; i >= 0; --i) {
        r <<= 1;
        r.lo |= (this->lo & (1ull << i)) > 0 ? 1 : 0;
        if (r >= other) {
          r -= other;
          q.lo |= 1ull << i;
        }
      }
      this->lo = r.lo;
      this->hi = r.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator+(uint128 a, uint128 b) { return create(a.value + b.value); }
#else
    friend inline uint128 operator+(uint128 a, uint128 b) { return a += b; }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator-(uint128 a, uint128 b) { return create(a.value - b.value); }
#else
    friend inline uint128 operator-(uint128 a, uint128 b) { return a -= b; }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator*(uint128 a, uint128 b) { return create(a.value * b.value); }
#else
    friend inline uint128 operator*(uint128 a, uint128 b) { return a *= b; }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator/(uint128 a, uint128 b) { return create(a.value / b.value); }
#else
    friend inline uint128 operator/(uint128 a, uint128 b) { return a /= b; }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator%(uint128 a, uint128 b) { return create(a.value % b.value); }
#else
    friend inline uint128 operator%(uint128 a, uint128 b) { return a %= b; }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator++() {
      ++value;
      return *this;
    }
#else
    uint128& operator++() {
      ++lo;
      if (lo == 0)
        ++hi;
      return *this;
    }
#endif

    uint128 operator++(int) {
      auto copy{*this};
      ++(*this);
      return copy;
    }

#if defined __SIZEOF_INT128__
    uint128& operator--() {
      --value;
      return *this;
    }
#else
    uint128& operator--() {
      if (lo == 0)
        --hi;
      --lo;
      return *this;
    }
#endif

    uint128 operator--(int) {
      auto copy{*this};
      --(*this);
      return copy;
    }

#if defined __SIZEOF_INT128__
    uint128& operator|=(uint128 other) {
      value |= other.value;
      return *this;
    }
#else
    uint128& operator|=(uint128 other) {
      lo |= other.lo;
      hi |= other.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator&=(uint128 other) {
      value &= other.value;
      return *this;
    }
#else
    uint128& operator&=(uint128 other) {
      lo &= other.lo;
      hi &= other.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator^=(uint128 other) {
      value ^= other.value;
      return *this;
    }
#else
    uint128& operator^=(uint128 other) {
      lo ^= other.lo;
      hi ^= other.hi;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator|(uint128 a, uint128 b) { return create(a.value | b.value); }
#else
    friend inline uint128 operator|(uint128 a, uint128 b) {
      return create(a.hi | b.hi, a.lo | b.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator&(uint128 a, uint128 b) { return create(a.value & b.value); }
#else
    friend inline uint128 operator&(uint128 a, uint128 b) {
      return create(a.hi & b.hi, a.lo & b.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator^(uint128 a, uint128 b) { return create(a.value ^ b.value); }
#else
    friend inline uint128 operator^(uint128 a, uint128 b) {
      return create(a.hi ^ b.hi, a.lo ^ b.lo);
    }
#endif

#if defined __SIZEOF_INT128__
    uint128 operator~() const { return create(~value); }
#else
    uint128 operator~() const { return create(~hi, ~lo); }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator>>=(int s) {
      value >>= s;
      return *this;
    }
#else
    uint128& operator>>=(int s) {
      if (s < 0)
        return (*this) <<= -s;
      uint128 result((this->hi >> s), (this->lo >> s) | (this->hi << (64 - s)));
      *this = result;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    uint128& operator<<=(int s) {
      value <<= s;
      return *this;
    }
#else
    uint128& operator<<=(int s) {
      if (s < 0)
        return (*this) >>= -s;
      uint128 result((this->hi << s) | (this->lo >> (64 - s)), this->lo << s);
      *this = result;
      return *this;
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator>>(uint128 a, int s) { return create(a.value >> s); }
#else
    friend inline uint128 operator>>(uint128 a, int s) {
      if (s < 0)
        return a << (-s);
      return uint128::create(a.hi >> s, (a.lo >> s) | (a.hi << (64 - s)));
    }
#endif

#if defined __SIZEOF_INT128__
    friend inline uint128 operator<<(uint128 a, int s) { return create(a.value << s); }
#else
    friend inline uint128 operator<<(uint128 a, int s) {
      if (s < 0)
        return a >> (-s);
      return uint128::create((a.hi << s) | (a.lo >> (64 - s)), a.lo << s);
    }
#endif
  };

}  // namespace trng

#endif
