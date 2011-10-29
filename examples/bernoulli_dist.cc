// Copyright (C) 2001-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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
#include <iomanip>
#include <vector>
#include <trng/lcg64.hpp>
#include <trng/bernoulli_dist.hpp>

typedef enum { head=0, tail=1 } coin;

int main() {
  // discrete distribution object
  trng::bernoulli_dist<coin> biased_coin(0.51, head, tail);
  // random number generator
  trng::lcg64 r;
  // draw some random numbers 
  std::vector<int> count(2, 0);
  const int samples=100000;
  for (int i=0; i<samples; ++i) {
    int x=biased_coin(r);  // draw a random number
    ++count[x];            // count
  }
  // print results
  std::cout << "value\t\tprobability\tcount\t\tempirical probability\n"
   	    << "=====\t\t===========\t=====\t\t=====================\n";
  for (std::vector<int>::size_type i=0; i<count.size(); ++i)
    std::cout << std::setprecision(3) 
	      << i << "\t\t"
	      << biased_coin.pdf(static_cast<coin>(i)) << "\t\t"
 	      << count[i] << "\t\t"
 	      << static_cast<double>(count[i])/samples << '\n';
  return EXIT_SUCCESS;
}
