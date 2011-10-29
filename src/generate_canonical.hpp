// Copyright (C) 2000-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#if !(defined TRNG_GENERATE_CANONICAL_HPP)

#define TRNG_GENERATE_CANONICAL_HPP

#include <trng/limits.hpp>
#include <trng/math.hpp>
#include <trng/utility.hpp>

namespace trng {
  
  template<typename result_type, typename R>
  result_type generate_canonical(R &);

  namespace detail {

    template<typename result_type, typename R, typename category>
    result_type generate_canonical_impl(R &, result_type, category);
    
    class integer_tag { };
    class float_tag { };
    
    template<bool>
    struct integer_float_traits;
    
    template<>
    struct integer_float_traits<false> {
      typedef float_tag cat;
    };
    
    template<>
    struct integer_float_traits<true> {
      typedef integer_tag cat;
    };
    
    template<typename result_type, typename R>
    inline result_type generate_canonical_impl(R &r, result_type, float_tag) {
      return utility::uniformoo<result_type>(r);
    }
  
    template<typename result_type, typename R>
    inline result_type generate_canonical_impl(R &r, result_type, integer_tag) {
      return static_cast<result_type>(math::floor(utility::uniformco<double>(r)*(static_cast<double>(R::max)-static_cast<double>(R::min)+1.0)));
    }
    
  }
  
  template<typename result_type, typename R>
  result_type generate_canonical(R &g) {
    return detail::generate_canonical_impl(g, result_type(), 
					   typename detail::integer_float_traits<math::numeric_limits<result_type>::is_integer>::cat());
  }
  
}

#endif
