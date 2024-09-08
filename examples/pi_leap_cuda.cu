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
#include <vector>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

__global__ void parallel_pi(long samples, trng::yarn2 *rx, trng::yarn2 *ry, long *in) {
  long rank = threadIdx.x;
  long size = blockDim.x;
  trng::uniform01_dist<float> u;  // random number distribution
  in[rank] = 0;                   // local number of points in circle
  for (long i = rank * samples / size; i < (rank + 1) * samples / size; ++i) {
    const float x = u(rx[rank]), y = u(ry[rank]);  // choose random x- and y-coordinates
    if (x * x + y * y <= 1)                        // is point in circle?
      ++in[rank];                                  // increase thread-local counter
  }
}

int main(int argc, char *argv[]) {
  const long samples{1000000l};            // total number of points in square
  const int size{128};                     // number of threads
  trng::yarn2 *rx{new trng::yarn2[size]};  // random number engines
  trng::yarn2 *ry{new trng::yarn2[size]};  // random number engines
  for (int rank{0}; rank < size; ++rank) {
    rx[rank].split(2, 0);        // choose sub-stream no. 0 out of 2 streams
    ry[rank].split(2, 1);        // choose sub-stream no. 1 out of 2 streams
    rx[rank].split(size, rank);  // choose sub-stream no. rank out of size streams
    ry[rank].split(size, rank);  // choose sub-stream no. rank out of size streams
  }
  // copy random number engines to CUDA device
  trng::yarn2 *rx_device, *ry_device;
  cudaMalloc(&rx_device, size * sizeof(*rx_device));
  cudaMalloc(&ry_device, size * sizeof(*ry_device));
  cudaMemcpy(rx_device, rx, size * sizeof(*rx), cudaMemcpyHostToDevice);
  cudaMemcpy(ry_device, ry, size * sizeof(*ry), cudaMemcpyHostToDevice);
  // memory for thread local results
  long *in_device;
  cudaMalloc(&in_device, size * sizeof(*in_device));
  // start parallel Monte Carlo
  parallel_pi<<<1, size>>>(samples, rx_device, ry_device, in_device);
  // gather results
  std::vector<long> in(size);
  cudaMemcpy(in.data(), in_device, size * sizeof(*in_device), cudaMemcpyDeviceToHost);
  cudaFree(rx_device);
  cudaFree(ry_device);
  long sum{0};
  for (int rank{0}; rank < size; ++rank)
    sum += in[rank];
  // print result
  std::cout << "pi = " << 4.0 * sum / samples << std::endl;
  return EXIT_SUCCESS;
}
