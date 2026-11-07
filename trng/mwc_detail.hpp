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

#if !(defined TRNG_MWC_DETAIL_HPP)

#define TRNG_MWC_DETAIL_HPP

#include <trng/cuda.hpp>
#include <trng/int_types.hpp>

namespace trng {

  // Internal helpers of the scrambled multiply-with-carry engines.  The state of such an
  // engine is an integer s modulo m = A * b^k - 1.  Each step multiplies s by b^-1 mod m,
  // such that a jump by n steps is a multiplication by (b^-1)^n mod m.  The moduli exceed
  // 64 bits, therefore integers are stored as arrays of 32-bit limbs (least significant
  // limb first) and the arithmetic needs neither division nor 128-bit integer types.
  namespace mwc_detail {

    // hi * 2^64 + lo = a * b + c
    TRNG_CUDA_ENABLE
    inline void mul_add(const uint64_t a, const uint64_t b, const uint64_t c, uint64_t &hi,
                        uint64_t &lo) {
#if defined __SIZEOF_INT128__ && !(defined __CUDACC__) && !(defined __HIPCC__)
      __extension__ typedef unsigned __int128 uint128_native;
      const uint128_native t{static_cast<uint128_native>(a) * b + c};
      lo = static_cast<uint64_t>(t);
      hi = static_cast<uint64_t>(t >> 64u);
#else
      const uint64_t mask{0xffffffffu};
      const uint64_t a0{a & mask}, a1{a >> 32u}, b0{b & mask}, b1{b >> 32u};
      const uint64_t p00{a0 * b0}, p01{a0 * b1}, p10{a1 * b0}, p11{a1 * b1};
      const uint64_t mid{(p00 >> 32u) + (p01 & mask) + (p10 & mask)};
      lo = (p00 & mask) | (mid << 32u);
      hi = p11 + (p01 >> 32u) + (p10 >> 32u) + (mid >> 32u);
      lo += c;
      if (lo < c)
        ++hi;
#endif
    }

    template<int N>
    struct limbs {
      uint32_t w[N];
    };

    // a += b, returns the carry
    template<int N>
    TRNG_CUDA_ENABLE uint32_t add(limbs<N> &a, const limbs<N> &b) {
      uint64_t carry{0};
      for (int i{0}; i < N; ++i) {
        const uint64_t t{static_cast<uint64_t>(a.w[i]) + b.w[i] + carry};
        a.w[i] = static_cast<uint32_t>(t);
        carry = t >> 32u;
      }
      return static_cast<uint32_t>(carry);
    }

    // a -= b, modulo 2^(32 N)
    template<int N>
    TRNG_CUDA_ENABLE void sub(limbs<N> &a, const limbs<N> &b) {
      uint64_t borrow{0};
      for (int i{0}; i < N; ++i) {
        const uint64_t t{static_cast<uint64_t>(a.w[i]) - b.w[i] - borrow};
        a.w[i] = static_cast<uint32_t>(t);
        borrow = (t >> 32u) & 1u;
      }
    }

    // a >= b ?
    template<int N>
    TRNG_CUDA_ENABLE bool geq(const limbs<N> &a, const limbs<N> &b) {
      for (int i{N - 1}; i >= 0; --i) {
        if (a.w[i] > b.w[i])
          return true;
        if (a.w[i] < b.w[i])
          return false;
      }
      return true;
    }

    // a = (a + b) mod m, requires a < m and b < m
    template<int N>
    TRNG_CUDA_ENABLE void add_mod(limbs<N> &a, const limbs<N> &b, const limbs<N> &m) {
      const uint32_t carry{add(a, b)};
      if (carry != 0 or geq(a, m))
        sub(a, m);
    }

    // a = (a * b) mod m, requires a < m and b < m (shift-and-add, most significant bit of
    // a first)
    template<int N>
    TRNG_CUDA_ENABLE void mul_mod(limbs<N> &a, const limbs<N> &b, const limbs<N> &m) {
      limbs<N> r;
      for (int i{0}; i < N; ++i)
        r.w[i] = 0;
      for (int i{N - 1}; i >= 0; --i)
        for (int k{31}; k >= 0; --k) {
          add_mod(r, r, m);
          if (((a.w[i] >> k) & 1u) != 0)
            add_mod(r, b, m);
        }
      a = r;
    }

    // r = base^n mod m
    template<int N>
    TRNG_CUDA_ENABLE void pow_mod(limbs<N> &r, const limbs<N> &base, unsigned long long n,
                                  const limbs<N> &m) {
      limbs<N> p{base};
      for (int i{0}; i < N; ++i)
        r.w[i] = 0;
      r.w[0] = 1;
      while (n != 0) {
        if ((n & 1u) != 0)
          mul_mod(r, p, m);
        n >>= 1u;
        if (n != 0)
          mul_mod(p, p, m);
      }
    }

    // r = base^(2^s) mod m
    template<int N>
    TRNG_CUDA_ENABLE void pow2_mod(limbs<N> &r, const limbs<N> &base, unsigned s,
                                   const limbs<N> &m) {
      r = base;
      for (unsigned i{0}; i < s; ++i)
        mul_mod(r, r, m);
    }

  }  // namespace mwc_detail

}  // namespace trng

#endif
