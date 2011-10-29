// Copyright (C) 2006 Heiko Bauke <heiko.bauke@physik.uni-magdeburg.de>
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
      result_type t=(1.0+static_cast<result_type>(r()-R::min))/
	(2.0+static_cast<result_type>(R::max-R::min));
      if (t==0.0)
	return math::numeric_limits<result_type>::epsilon();
      if (t==1.0)
	return 1.0-math::numeric_limits<result_type>::epsilon();
      return t;
    }
  
    template<typename result_type, typename R>
    inline result_type generate_canonical_impl(R &r, result_type, integer_tag) {
      double t=(static_cast<double>(r()-R::min))/
	(1.0+static_cast<double>(R::max-R::min))*
	(static_cast<double>(math::numeric_limits<result_type>::max())-
	 static_cast<double>(math::numeric_limits<result_type>::min()))+
	static_cast<double>(math::numeric_limits<result_type>::min());
      return static_cast<result_type>(t);
    }
    
  }
  
  template<typename result_type, typename R>
  result_type generate_canonical(R &g) {
    return detail::generate_canonical_impl(g, result_type(), 
					   typename detail::integer_float_traits<math::numeric_limits<result_type>::is_integer>::cat());
  }
  
}

#endif
