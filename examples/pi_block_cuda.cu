// Copyright (c) 2000-2022, Heiko Bauke
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
#include <trng/yarn5s.hpp>
#include <trng/uniform01_dist.hpp>

__global__ void parallel_pi(long samples, long *in, trng::yarn5s r) {
  long rank = threadIdx.x;
  long size = blockDim.x;
  r.jump(2 * (rank * samples / size));  // jump ahead
  trng::uniform01_dist<float> u;        // random number distribution
  in[rank] = 0;                         // local number of points in circle
  for (long i = rank * samples / size; i < (rank + 1) * samples / size; ++i) {
    const float x = u(r), y = u(r);  // choose random x- and y-coordinates
    if (x * x + y * y <= 1)          // is point in circle?
      ++in[rank];                    // increase thread-local counter
  }
}

int main(int argc, char *argv[]) {
  const long samples{1000000l};  // total number of points in square
  const int size{128};           // number of threads
  long *in_device;
  cudaMalloc(&in_device, size * sizeof(*in_device));
  trng::yarn5s r;
  // start parallel Monte Carlo
  parallel_pi<<<1, size>>>(samples, in_device, r);
  // gather results
  std::vector<long> in(size);
  cudaMemcpy(in.data(), in_device, size * sizeof(*in_device), cudaMemcpyDeviceToHost);
  cudaFree(in_device);
  long sum{0};
  for (int rank{0}; rank < size; ++rank)
    sum += in[rank];
  // print result
  std::cout << "pi = " << 4.0 * sum / samples << std::endl;
  return EXIT_SUCCESS;
}
