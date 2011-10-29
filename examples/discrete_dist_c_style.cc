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
#include <iomanip>
#include <trng/lcg64.hpp>
#include <trng/discrete_dist.hpp>

int main() {
  const int n=6;
  // stores relative probabilities
  double p[6]={ 1, 3.25, 5, 6.5, 7, 2 };
  // discrete distribution object
  trng::discrete_dist dist(p, p+n);
  // random number generator
  trng::lcg64 r;
  // draw some random numbers 
  int count[n]={ 0, 0, 0, 0, 0, 0 };
  const int samples=10000;
  for (int i=0; i<samples; ++i) {
    int x=dist(r);  // draw a random number
    ++count[x];     // count
  }
  // print results
  std::cout << "value\t\tprobability\tcount\t\tempirical probability\n"
	    << "=====\t\t===========\t=====\t\t=====================\n";
  for (int i=0; i<n; ++i) {
    std::cout << std::setprecision(3) 
	      << i << "\t\t"
	      << dist.pdf(i) << "\t\t"
	      << count[i] << "\t\t"
	      << static_cast<double>(count[i])/samples << '\n';
  }
  return EXIT_SUCCESS;
}
