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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <trng/lcg64.hpp>
#include <trng/mrg2.hpp>
#include <trng/mrg3.hpp>
#include <trng/mrg3s.hpp>
#include <trng/mrg4.hpp>
#include <trng/mrg5.hpp>
#include <trng/mrg5s.hpp>
#include <trng/yarn2.hpp>
#include <trng/yarn3.hpp>
#include <trng/yarn3s.hpp>
#include <trng/yarn4.hpp>
#include <trng/yarn5.hpp>
#include <trng/yarn5s.hpp>


template<typename rng_type>
void implemtation_test() {
  bool test_passed;

  std::cout << "Testing generator " << rng_type::name() << "\n"
	    << "============================================\n\n";
  // NumberGenerator concept
  {
    test_passed=true;
    if (!std::numeric_limits<typename rng_type::result_type>::is_specialized)
      test_passed=false;
    rng_type R;
    typename rng_type::result_type t=R();
    std::cout << "NumberGenerator concept:\t\t" << (test_passed ? "ok" : "fail") << "\n";
  }
  // UniformRandomNumberGenerator
  {
    test_passed=true;
    rng_type R;
    typename rng_type::result_type min=R.min, max=R.max;
    if (min>max)
      test_passed=false;
    std::cout << "UniformRandomNumberGenerator concept:\t" << (test_passed ? "ok" : "fail") << "\n";
    std::cout << "result_type is " 
	      << (std::numeric_limits<typename rng_type::result_type>::is_integer ? "" : "not")
	      << "an integer type\n";
  }
  // PseudoRandomNumberGenerator concept
  {
    test_passed=true;
    rng_type R_init;
    rng_type R1(R_init), R2, R3;
    R2=R1;
    if (R1()!=R2())
      test_passed=false;
    R1.seed(R_init);
    R2=R1;
    if (R1()!=R2())
      test_passed=false;
    R1.seed();
    if (R1()!=R3())
      test_passed=false;
    R1.seed(R_init);
    std::stringstream str;
    str << R3;
    str >> R1;
    if (R1()!=R3())
      test_passed=false;
    std::cout << "PseudoRandomNumberGenerator concept:\t" << (test_passed ? "ok" : "fail") << "\n";
  }
  
  std::cout << "\n\n";
}


int main() {
  implemtation_test<trng::lcg64>();  
  implemtation_test<trng::mrg2>();  
  implemtation_test<trng::mrg3>();  
  implemtation_test<trng::mrg3s>();  
  implemtation_test<trng::mrg4>();  
  implemtation_test<trng::mrg5>();  
  implemtation_test<trng::mrg5s>();  
  implemtation_test<trng::yarn2>();  
  implemtation_test<trng::yarn3>();  
  implemtation_test<trng::yarn3s>();  
  implemtation_test<trng::yarn4>();  
  implemtation_test<trng::yarn5>();  
  implemtation_test<trng::yarn5s>();  
  return EXIT_SUCCESS;
}
