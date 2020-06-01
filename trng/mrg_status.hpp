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

#if !(defined TRNG_MRG_STATUS_HPP)
#define TRNG_MRG_STATUS_HPP

#include <ostream>
#include <istream>
#include <algorithm>

namespace trng {

  template<typename result_type, int n, typename F>
  class mrg_status {
  protected:
    result_type r[n];

  public:
    mrg_status() {
      r[0] = 0;
      for (int i{1}; i < n; ++i)
        r[i] = 1;
    }
    explicit mrg_status(result_type r1, result_type r2) : r{r1, r2} {
      static_assert(n == 2, "dimension must be 2");
    }
    explicit mrg_status(result_type r1, result_type r2, result_type r3) : r{r1, r2, r3} {
      static_assert(n == 3, "dimension must be 3");
    }
    explicit mrg_status(result_type r1, result_type r2, result_type r3, result_type r4)
        : r{r1, r2, r3, r4} {
      static_assert(n == 4, "dimension must be 4");
    }
    explicit mrg_status(result_type r1, result_type r2, result_type r3, result_type r4,
                        result_type r5)
        : r{r1, r2, r3, r4, r5} {
      static_assert(n == 5, "dimension must be 5");
    }

    friend F;

    // Equality comparable concept
    friend bool operator==(const mrg_status &S1, const mrg_status &S2) {
      return std::equal(S1.r, S1.r + n, S2.r);
    }

    friend bool operator!=(const mrg_status &P1, const mrg_status &P2) { return not(P1 == P2); }


    // Streamable concept
    template<typename char_t, typename traits_t>
    friend std::basic_ostream<char_t, traits_t> &operator<<(
        std::basic_ostream<char_t, traits_t> &out, const mrg_status &S) {
      std::ios_base::fmtflags flags(out.flags());
      out.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      out << '(';
      for (int i{0}; i < n; ++i) {
        out << S.r[i];
        if (i + 1 < n)
          out << ' ';
      }
      out << ')';
      out.flags(flags);
      return out;
    }

    template<typename char_t, typename traits_t>
    friend std::basic_istream<char_t, traits_t> &operator>>(
        std::basic_istream<char_t, traits_t> &in, mrg_status &S) {
      mrg_status S_new;
      std::ios_base::fmtflags flags(in.flags());
      in.flags(std::ios_base::dec | std::ios_base::fixed | std::ios_base::left);
      in >> utility::delim('(');
      for (int i{0}; i < n; ++i) {
        in >> S_new.r[i];
        if (i + 1 < n)
          in >> utility::delim(' ');
      }
      in >> utility::delim(')');
      if (in)
        S = S_new;
      in.flags(flags);
      return in;
    }
  };

}  // namespace trng

#endif
