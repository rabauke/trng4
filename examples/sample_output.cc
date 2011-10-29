// Copyright (C) 2001-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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
#include <exception>
#include <trng/config.hpp>
#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
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
// #include <trng/mt19937.hpp>
// #include <trng/mt19937_64.hpp>
#include <trng/lagfib2xor.hpp>
#include <trng/lagfib2plus.hpp>
#include <trng/lagfib4xor.hpp>
#include <trng/lagfib4plus.hpp>


template<typename R>
void sample_output(R &r, std::string name) {
  while (name.length()<32)
    name+=' ';
  std::cout << name;
  for (int i=0; i<15; ++i)
    std::cout << r() << '\t';
  std::cout << r() << '\n';
}


int main() {
  try {
    { trng::lcg64       r;  sample_output(r, "trng::lcg64"); }
    { trng::lcg64_shift r;  sample_output(r, "trng::lcg64_shift"); }
    { trng::mrg2        r;  sample_output(r, "trng::mrg2"); }
    { trng::mrg3        r;  sample_output(r, "trng::mrg3"); }
    { trng::mrg3s       r;  sample_output(r, "trng::mrg3s"); }
    { trng::mrg4        r;  sample_output(r, "trng::mrg4"); }
    { trng::mrg5        r;  sample_output(r, "trng::mrg5"); }
    { trng::mrg5s       r;  sample_output(r, "trng::mrg5s"); }
    { trng::yarn2       r;  sample_output(r, "trng::yarn2"); }
    { trng::yarn3       r;  sample_output(r, "trng::yarn3"); }
    { trng::yarn3s      r;  sample_output(r, "trng::yarn3s"); }
    { trng::yarn4       r;  sample_output(r, "trng::yarn4"); }
    { trng::yarn5       r;  sample_output(r, "trng::yarn5"); }
    { trng::yarn5s      r;  sample_output(r, "trng::yarn5s"); }
    // { trng::mt19937     r;  sample_output(r, "trng::mt19937"); }
    // { trng::mt19937_64  r;  sample_output(r, "trng::mt19937_64"); }
    { trng::lagfib2xor_19937_ull  r;  sample_output(r, "trng::lagfib2xor_19937_ull"); }
    { trng::lagfib4xor_19937_ull  r;  sample_output(r, "trng::lagfib4xor_19937_ull"); }
    { trng::lagfib4plus_19937_ull r;  sample_output(r, "trng::lagfib2plus_19937_ull"); }
    { trng::lagfib2plus_19937_ull r;  sample_output(r, "trng::lagfib4plus_19937_ull"); }
  } 
  catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
  }
  catch (...) {
    std::cerr << "something went wrong" << std::endl;
  }
  return EXIT_SUCCESS;
}
