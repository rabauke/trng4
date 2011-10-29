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
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

int main(int argc, char *argv[]) {
  const long samples=1000000l;          // total number of points in square
  long in=0l;                           // no points in circle
  trng::yarn2 r;                        // random number engine
  trng::uniform01_dist u;               // random number distribution
  // throw random points into square 
  for (long i=0; i<samples; ++i) {
    double x=u(r), y=u(r);              // choose random x- and y-coordinates
    if (x*x+y*y<=1.0)                   // is point in circle?
      ++in;                             // increase counter
  }
  std::cout << "pi = " << 4.0*in/samples << std::endl;
  return EXIT_SUCCESS;
}
