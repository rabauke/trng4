// Copyright (c) 2000-2014, Heiko Bauke
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
//  

#if !(defined TRNG_UNIFORMXX_HPP)
#define TRNG_UNIFORMXX_HPP

#include <cstddef>
#include <cfloat>
#include <trng/config.hpp>
#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/static_assertion.hpp>

namespace trng {

  namespace utility {

    template<unsigned long long top, unsigned int count=0>
    struct Bits {
      enum { result=Bits<top/2, count+1>::result }; 
    };

    template<unsigned int count>
    struct Bits<0ULL, count> {
      enum { result=count };
    };

    //------------------------------------------------------------------

    template<unsigned long long top, unsigned int count=0>
    struct Holes {
      enum { result=Holes<top/2, count + (~top & 1)>::result };
    };

    template<unsigned int count>
    struct Holes<0ULL, count> {
      enum { result=count };
    };

    //------------------------------------------------------------------

    template<typename F>
    struct epsilon;

    template<>
    struct epsilon<float> {
      TRNG_CUDA_ENABLE
#if !(defined __CUDA_ARCH__)
      static 
#endif
      float val() {
	return FLT_EPSILON;
      }
    };

    template<>
    struct epsilon<double> {
      TRNG_CUDA_ENABLE
#if !(defined __CUDA_ARCH__)
      static 
#endif
      double val() {
	return DBL_EPSILON;
      }
    };
    
    template<>
    struct epsilon<long double> {
      TRNG_CUDA_ENABLE
#if !(defined __CUDA_ARCH__)
      static 
#endif
      long double val() {
	return LDBL_EPSILON;
      }
    };    

    //------------------------------------------------------------------

    // With basic optimizations enabled, modern C++ compilers can reduce
    // all the public routines herein down to small inline code sequences.
    // They should also collapse the size (sizeof(u01xx_traits<...>) to 1.
    template<typename return_type, std::size_t requested_bits, typename prng_t>
    class u01xx_traits {
      typedef return_type ret_t;
      typedef typename prng_t::result_type result_type;

      // Casting up from "simpler" types may yield better low level code sequences.
#if __cplusplus>=201103L
      static const result_type domain_max0=prng_t::max()-prng_t::min();
#else
      static const result_type domain_max0=prng_t::max-prng_t::min;
#endif
      static const unsigned int domain_bits=Bits<domain_max0>::result;
      static const unsigned int domain_full_bits=domain_bits-(Holes<domain_max0>::result > 0);
      static const bool int_ok=domain_bits < static_cast<unsigned int>(math::numeric_limits<unsigned int>::digits);
      static const bool long_ok=domain_bits < static_cast<unsigned int>(math::numeric_limits<unsigned long>::digits);
      static const bool long_long_ok=domain_bits < static_cast<unsigned int>(math::numeric_limits<unsigned long long>::digits);

      static const unsigned int ret_bits=math::numeric_limits<ret_t>::digits;
      static const unsigned int bits0=(requested_bits < ret_bits) ? requested_bits : ret_bits;
      static const unsigned int bits=(bits0 < 1) ? 1 : bits0;
      static const std::size_t calls_needed=(bits/domain_full_bits) + ((bits % domain_full_bits) != 0);

      // a ((long long)(val >> 1)) cast may give us better performance (~4X using lcg on Core2 x86-64)
      // the lost bit will not usually be missed as it is ~30dB down typically (53 bit mantissa from
      // a 63 rather than 64 bit integer variate)
      static const bool use_ll_of_shifted=not long_long_ok and (domain_max0 >> 1)==(~0ULL >> 1) and bits<domain_bits;
      static const result_type domain_max=use_ll_of_shifted ? (domain_max0 >> 1) : domain_max0;

      TRNG_CUDA_ENABLE
      static ret_t addin(prng_t &r) {
#if __cplusplus>=201103L
	result_type x=r()-prng_t::min();
#else
	result_type x=r()-prng_t::min;
#endif
	if (int_ok) 
	  return static_cast<int>(x);
	else if (long_ok) 
	  return static_cast<long>(x);
	else if (long_long_ok) 
	  return static_cast<long long>(x);
	else if (use_ll_of_shifted) 
	  return static_cast<long long>(x>>1);
	else 
	  return x;
      }
      TRNG_CUDA_ENABLE
      static ret_t variate_max() {
	const ret_t scale_per_step(ret_t(domain_max) + 1);
	const ret_t max_addin(domain_max);
	ret_t ret(max_addin);
	for (std::size_t i(1); i<calls_needed; ++i)
	  ret=ret*scale_per_step+max_addin;
	return ret;
      }
      TRNG_CUDA_ENABLE
      static ret_t variate(prng_t &r) {
#if !(defined __CUDA_ARCH__)
#if __cplusplus>=201103L
	detail::static_assertion<(prng_t::min()>=0 and prng_t::max()>prng_t::min())>(); // min and/or max out of spec?
	detail::static_assertion<(prng_t::max()-prng_t::min()) <= ~0ULL>(); // Bits, Holes incorrect otherwise
#else
	detail::static_assertion<(prng_t::min>=0 and prng_t::max>prng_t::min)>(); // min and/or max out of spec?
	detail::static_assertion<(prng_t::max-prng_t::min) <= ~0ULL>(); // Bits, Holes incorrect otherwise
#endif
	detail::static_assertion<not math::numeric_limits<return_type>::is_integer>();
	detail::static_assertion<(calls_needed>0 and calls_needed<=bits)>();
#endif
	const ret_t scale_per_step(ret_t(domain_max) + 1);
	ret_t ret(addin(r));
	for (std::size_t i(1); i<calls_needed; ++i)
	  ret=ret*scale_per_step+addin(r);
	return ret;
      }
      TRNG_CUDA_ENABLE
      static ret_t eps() {
#if defined __CUDA_ARCH__
	epsilon<ret_t> EPS;
        const ret_t native_eps(EPS.val());
#else
        const ret_t native_eps(epsilon<ret_t>::val());
#endif
        const ret_t domain_eps(ret_t(1)/domain_max);
	return native_eps>=domain_eps or requested_bits!=1 ? native_eps : domain_eps;
      }
      TRNG_CUDA_ENABLE
      static ret_t cc_norm() {
	return ret_t(1)/variate_max();
      }
      TRNG_CUDA_ENABLE
      static ret_t co_norm() {
	return cc_norm()*(ret_t(1)-eps());
      }
      TRNG_CUDA_ENABLE
      static ret_t oo_norm() {
	return cc_norm()*(ret_t(1)-2*eps());
      }
    public:
      TRNG_CUDA_ENABLE
      static return_type cc(prng_t &r) {
        const bool division_required(variate_max()*cc_norm()!=1);
        return division_required ? variate(r)/variate_max() : variate(r)*cc_norm();
      }
      TRNG_CUDA_ENABLE
      static return_type co(prng_t &r) { 
	return variate(r)*co_norm(); 
      }
      TRNG_CUDA_ENABLE
      static return_type oc(prng_t &r) { 
	return ret_t(1)-co(r);	
      }
      TRNG_CUDA_ENABLE
      static return_type oo(prng_t &r) { 
	return variate(r)*oo_norm()+eps(); 
      }
    };

    template<typename ReturnType, std::size_t bits, typename UniformRandomNumberGenerator>
    TRNG_CUDA_ENABLE
    ReturnType generate_canonical(UniformRandomNumberGenerator &u) {
      return u01xx_traits<ReturnType, bits, UniformRandomNumberGenerator>::co(u);
    }

    template<typename ReturnType, typename PrngType>
    TRNG_CUDA_ENABLE
    inline ReturnType uniformcc(PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::cc(r); 
    }

    template<typename ReturnType, typename PrngType>
    TRNG_CUDA_ENABLE
    inline ReturnType uniformco(PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::co(r); 
    }

    template<typename ReturnType, typename PrngType>
    TRNG_CUDA_ENABLE
    inline ReturnType uniformoc(PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::oc(r); 
    }

    template<typename ReturnType, typename PrngType>
    TRNG_CUDA_ENABLE
    inline ReturnType uniformoo(PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::oo(r); 
    }

  }

}

#endif
