// Copyright (c) 2000-2024, Heiko Bauke
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.
//
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.

#include <cstdlib>
#include <iostream>
#include "mpi.h"
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

int main(int argc, char *argv[]) {
  const long samples{1000000l};  // total number of points in square
  MPI_Init(&argc, &argv);        // initialise MPI environment
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);  // get total number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get rank of current process
  long in{0};                            // number of points in circle
  trng::yarn2 r;                         // random number engine
  trng::uniform01_dist<> u;              // random number distribution
  r.jump(2 * (rank * samples / size));   // jump ahead
  // throw random points into square and distribute workload over all processes
  for (long i{rank * samples / size}; i < (rank + 1) * samples / size; ++i) {
    const double x{u(r)}, y{u(r)};  // choose random x- and y-coordinates
    if (x * x + y * y <= 1.0)       // is point in circle?
      ++in;                         // increase counter
  }
  // calculate sum of all local variables 'in' and storre result in 'in_all' on process 0
  long in_all;
  MPI_Reduce(&in, &in_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0)  // print result
    std::cout << "pi = " << 4.0 * in_all / samples << std::endl;
  MPI_Finalize();  // quit MPI
  return EXIT_SUCCESS;
}
