// Copyright (c) 2000-2018, Heiko Bauke
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

#if !(defined TRNG_INT_TYPES_HPP)

#define TRNG_INT_TYPES_HPP

#if __cplusplus>=201103L

#include<cstdint>

namespace trng {

  typedef std::int8_t int8_t;
  typedef std::int16_t int16_t;
  typedef std::int32_t int32_t;
  typedef std::int64_t int64_t;

  typedef std::uint8_t uint8_t;
  typedef std::uint16_t uint16_t;
  typedef std::uint32_t uint32_t;
  typedef std::uint64_t uint64_t;

}

#else

#include<climits>

#if CHAR_BIT!=8
#error "unsupported platform, type char must be of exactly 8 bits"
#endif

namespace trng {

  namespace detail {

    class invalid_type;
  
    template<bool, typename T1, typename T2=invalid_type>
    class if_then_else;
  
    template<typename T1, typename T2>
    class if_then_else<true, T1, T2> {
    public:
      typedef T1 type;
    };
  
    template<typename T1, typename T2>
    class if_then_else<false, T1, T2> {
    public:
      typedef T2 type;
    };
  
    template<int bits>
    class signed_integer {
    public:
      typedef 
      typename if_then_else<(sizeof(signed char)*CHAR_BIT>=bits), signed char, 
	typename if_then_else<(sizeof(signed short int)*CHAR_BIT>=bits), signed short int, 
	typename if_then_else<(sizeof(signed int)*CHAR_BIT>=bits), signed int, 
	typename if_then_else<(sizeof(signed long int)*CHAR_BIT>=bits), signed long int, 
	typename if_then_else<(sizeof(signed long long int)*CHAR_BIT>=bits), signed long long int>::type 
	>::type 
	>::type 
	>::type 
	>::type type;
    };

    template<int bits>
    class unsigned_integer {
    public:
      typedef 
      typename if_then_else<(sizeof(unsigned char)*CHAR_BIT>=bits), unsigned char, 
	typename if_then_else<(sizeof(unsigned short int)*CHAR_BIT>=bits), unsigned short int, 
	typename if_then_else<(sizeof(unsigned int)*CHAR_BIT>=bits), unsigned int, 
	typename if_then_else<(sizeof(unsigned long int)*CHAR_BIT>=bits), unsigned long int, 
	typename if_then_else<(sizeof(unsigned long long int)*CHAR_BIT>=bits), unsigned long long int>::type 
	>::type 
	>::type 
	>::type 
	>::type type;
    };

  }

  typedef detail::signed_integer<8>::type int8_t;
  typedef detail::signed_integer<16>::type int16_t;
  typedef detail::signed_integer<32>::type int32_t;
  typedef detail::signed_integer<64>::type int64_t;

  typedef detail::unsigned_integer<8>::type uint8_t;
  typedef detail::unsigned_integer<16>::type uint16_t;
  typedef detail::unsigned_integer<32>::type uint32_t;
  typedef detail::unsigned_integer<64>::type uint64_t;

}

#endif

#endif
