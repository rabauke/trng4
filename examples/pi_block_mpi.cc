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

#include <trng/config.hpp>
#if defined TRNG_HAVE_MPI

#include <cstdlib>
#include <iostream>
#include "mpi.h"
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

int main(int argc, char *argv[]) {
  const long samples=1000000l;          // total number of points in square
  trng::yarn2 r;                        // random number engine
  MPI::Init(argc, argv);                // initialise MPI environment
  int size=MPI::COMM_WORLD.Get_size();  // get total number of processes
  int rank=MPI::COMM_WORLD.Get_rank();  // get rank of current process
  long in=0l;                           // number of points in circle
  trng::uniform01_dist<> u;             // random number distribution
  r.jump(2*(rank*samples/size));        // jump ahead
  // throw random points into square and distribute workload over all processes
  for (long i=rank*samples/size; i<(rank+1)*samples/size; ++i) {
    double x=u(r), y=u(r);              // choose random x- and y-coordinates
    if (x*x+y*y<=1.0)                   // is point in circle?
      ++in;                             // increase counter
  }
  // calculate sum of all local variables 'in' and storre result in 'in_all' on process 0
  long in_all;
  MPI::COMM_WORLD.Reduce(&in, &in_all, 1, MPI::LONG, MPI::SUM, 0);
  if (rank==0)                          // print result 
    std::cout << "pi = " << 4.0*in_all/samples << std::endl;
  MPI::Finalize();                      // quit MPI
  return EXIT_SUCCESS;
}

#else

#include <cstdlib>
#include <iostream>

int main() {
  std::cerr << "Sorry, no MPI library installed on your computer.\n";
  return EXIT_FAILURE;
}

#endif
