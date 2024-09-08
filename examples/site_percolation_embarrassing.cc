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
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include "mpi.h"

const int number_of_realizations{1000};
const int Nx{250}, Ny{200};  // grid size
const int number_of_PRNs_per_sweep{Nx * Ny};
int site[Nx][Ny];      // lattice
const double P{0.46};  // occupation probability

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);  // initialize MPI environment
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);  // get total number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get rank of current process
  trng::yarn2 R;                         // random number engine
  trng::uniform01_dist<> u;              // random number distribution
  // skip random numbers that are consumed by other processes
  R.jump(rank * number_of_PRNs_per_sweep);
  for (int i{rank}; i < number_of_realizations; i += size) {
    // consume Nx * Ny pseudo-random numbers
    for (int x{0}; x < Nx; ++x)
      for (int y{0}; y < Ny; ++y)
        if (u(R) < P)
          site[x][y] = 1;  // site is occupied
        else
          site[x][y] = 0;  // site is not occupied
    // skip random numbers that are consumed by other processes
    R.jump((size - 1) * number_of_PRNs_per_sweep);
    // analyze lattice
    // ... source omitted
  }
  MPI_Finalize();  // quit MPI
  return EXIT_SUCCESS;
}
