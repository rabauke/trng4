// Copyright (c) 2000-2026, Heiko Bauke
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
#include <new>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include "mpi.h"

const int number_of_realizations{1000};
const int Nx{250}, Ny{200};  // grid size
const double P{0.46};        // occupation probability

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);  // initialize MPI environment
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);  // get total number of processes
  // create a two-dimensional Cartesian communicator
  int dims[2]{0, 0};             // number of processes in each domension
  int coords[2];                 // coordinates of current process within the grid
  int periods[2]{false, false};  // no periodic boundary conditions
  // calculate a balanced grid partitioning such that  size = dims[0] * dims[1]
  MPI_Dims_create(size, 2, dims);
  MPI_Comm Comm;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &Comm);
  int rank;
  MPI_Comm_rank(Comm, &rank);              // get rank of current process
  MPI_Cart_coords(Comm, rank, 2, coords);  // get coordinates of current process
  // determine section of current process
  int x0{coords[0] * Nx / dims[0]}, x1{(coords[0] + 1) * Nx / dims[0]}, Nxl{x1 - x0},
      y0{coords[1] * Ny / dims[1]}, y1{(coords[1] + 1) * Ny / dims[1]}, Nyl{y1 - y0};
  int *site{new int[Nxl * Nyl]};  // allocate memory to storre a sublattice
  trng::yarn2 R;                  // random number engine
  trng::uniform01_dist<> u;       // random number distribution
  // skip random numbers that are consumed by other processes
  R.jump(Nx * y0 + x0);
  for (int i{0}; i < number_of_realizations; ++i) {
    // consume Nxl * Nyl pseudo-random numbers
    int *s{site};
    for (int y{y0}; y < y1; ++y) {
      for (int x{x0}; x < x1; ++x) {
        if (u(R) < P)
          *s = 1;  // site is occupied
        else
          *s = 0;  // site is not occupied
        ++s;
      }
      // skip random numbers that are consumed by other processes
      R.jump(Nx - Nxl);
    }
    // skip random numbers that are consumed by other processes
    R.jump(Nx * (Ny - Nyl));
    // analyze lattice
    // ... source omitted
  }
  delete[] site;
  MPI_Finalize();  // quit MPI
  return EXIT_SUCCESS;
}
