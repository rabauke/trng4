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
// include TRNG header files
#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>

int main() {
  // random number engine
  trng::yarn2 R;
  // normal distribution with mean 6 and standard deviation 2
  trng::normal_dist<> normal(6.0, 2.0);
  // generate 1000 normal distributed random numbers
  for (int i=0; i<1000; ++i)
    std::cout << normal(R) << '\n';
  return EXIT_SUCCESS;
}
