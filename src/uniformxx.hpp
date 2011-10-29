// Copyright (C) 2000-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
// and Bruce Carneal
//  
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License in
// version 2 as published by the Free Software Foundation.
//  
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//  
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
//  

#if !(defined TRNG_UNIFORMXX_HPP)
#define TRNG_UNIFORMXX_HPP

#include <cstddef>
#include <trng/config.hpp>
#include <trng/limits.hpp>

#if TRNG_HAVE_BOOST
#include <boost/static_assert.hpp>
#else 
#define BOOST_STATIC_ASSERT(foo) 
#endif

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

    // With basic optimizations enabled, modern C++ compilers can reduce
    // all the public routines herein down to small inline code sequences.
    // They should also collapse the size (sizeof(u01xx_traits<...>) to 1.
    template<typename return_type, std::size_t requested_bits, typename prng_t>
    class u01xx_traits {
      typedef return_type ret_t;
      typedef typename prng_t::result_type result_type;

      BOOST_STATIC_ASSERT(prng_t::min >= 0 and prng_t::max > prng_t::min); // min and/or max out of spec?
      BOOST_STATIC_ASSERT((prng_t::max-prng_t::min) <= ~0ULL); // Bits, Holes incorrect otherwise
      BOOST_STATIC_ASSERT(!math::numeric_limits<return_type>::is_integer);

      // Casting up from "simpler" types may yield better low level code sequences.
      static const result_type domain_max0 = prng_t::max - prng_t::min;
      static const unsigned int domain_bits = Bits<domain_max0>::result;
      static const unsigned int domain_full_bits = domain_bits - (Holes<domain_max0>::result > 0);
      static const bool int_ok = domain_bits < math::numeric_limits<unsigned int>::digits;
      static const bool long_ok = domain_bits < math::numeric_limits<unsigned long>::digits;
      static const bool long_long_ok = domain_bits < math::numeric_limits<unsigned long long>::digits;

      static const unsigned int ret_bits = math::numeric_limits<ret_t>::digits;
      static const unsigned int bits0 = (requested_bits < ret_bits) ? requested_bits : ret_bits;
      static const unsigned int bits = (bits0 < 1) ? 1 : bits0;
      static const std::size_t calls_needed = (bits/domain_full_bits) + ((bits % domain_full_bits) != 0);
      BOOST_STATIC_ASSERT(calls_needed > 0 and calls_needed <= bits);

      // a ((long long)(val >> 1)) cast may give us better performance (~4X using lcg on Core2 x86-64)
      // the lost bit will not usually be missed as it is ~30dB down typically (53 bit mantissa from
      // a 63 rather than 64 bit integer variate)
      static const bool use_ll_of_shifted = not long_long_ok and (domain_max0 >> 1)==(~0ULL >> 1) and bits<domain_bits;
      static const result_type domain_max = use_ll_of_shifted ? (domain_max0 >> 1) : domain_max0;

      static ret_t addin(const prng_t &r) {
	result_type x=r()-prng_t::min;
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
      static ret_t variate_max() {
	const ret_t scale_per_step(ret_t(domain_max) + 1);
	const ret_t max_addin(domain_max);
	ret_t ret(max_addin);
	for (std::size_t i(1); i<calls_needed; ++i)
	  ret=ret*scale_per_step+max_addin;
	return ret;
      }
      static ret_t variate(const prng_t &r) {
	const ret_t scale_per_step(ret_t(domain_max) + 1);
	ret_t ret(addin(r));
	for (std::size_t i(1); i<calls_needed; ++i)
	  ret=ret*scale_per_step+addin(r);
	return ret;
      }
      static ret_t eps() {
        const ret_t native_eps(math::numeric_limits<ret_t>::epsilon());
        const ret_t domain_eps(ret_t(1)/domain_max);
	return native_eps>=domain_eps or requested_bits!=1 ? native_eps : domain_eps;
      }
      static ret_t cc_norm() {
	return ret_t(1)/variate_max();
      }
      static ret_t co_norm() {
	return cc_norm()*(ret_t(1)-eps());
      }
      static ret_t oo_norm() {
	return cc_norm()*(ret_t(1)-2*eps());
      }
    public:
      static return_type cc(const prng_t &r) {
        const bool division_required(variate_max()*cc_norm()!=1);
        return division_required ? variate(r)/variate_max() : variate(r)*cc_norm();
      }
      static return_type co(const prng_t &r) { 
	return variate(r)*co_norm(); 
      }
      static return_type oc(const prng_t &r) { 
	return ret_t(1)-co(r);	
      }
      static return_type oo(const prng_t &r) { 
	return variate(r)*oo_norm()+eps(); 
      }
    };

    template<typename ReturnType, std::size_t bits, typename UniformRandomNumberGenerator>
    ReturnType generate_canonical(UniformRandomNumberGenerator &u) {
      return u01xx_traits<ReturnType, bits, UniformRandomNumberGenerator>::co(u);
    }

    template<typename ReturnType, typename PrngType>
    inline ReturnType uniformcc(const PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::cc(r); 
    }

    template<typename ReturnType, typename PrngType>
    inline ReturnType uniformco(const PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::co(r); 
    }

    template<typename ReturnType, typename PrngType>
    inline ReturnType uniformoc(const PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::oc(r); 
    }

    template<typename ReturnType, typename PrngType>
    inline ReturnType uniformoo(const PrngType &r) { 
      return u01xx_traits<ReturnType, 1, PrngType>::oo(r); 
    }

  }

}

#endif
