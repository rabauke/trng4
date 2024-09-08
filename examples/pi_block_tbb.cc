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
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

class parallel_pi {
  trng::uniform01_dist<> u;  // random number distribution
  const trng::yarn2 &r;
  long in;

public:
  void operator()(const tbb::blocked_range<long> &range) {
    trng::yarn2 r_local(r);           // local copy of random number engine
    r_local.jump(2 * range.begin());  // jump ahead
    for (long i{range.begin()}; i != range.end(); ++i) {
      const double x{u(r_local)}, y{u(r_local)};  // choose random x- and y-coordinates
      if (x * x + y * y <= 1.0)                   // is point in circle?
        ++in;                                     // increase thread-local counter
    }
  }
  // join threds and counters
  void join(const parallel_pi &other) { in += other.in; }
  long in_circle() const { return in; }
  explicit parallel_pi(const trng::yarn2 &r) : r{r}, in{0} {}
  explicit parallel_pi(const parallel_pi &other, tbb::split) : r{other.r}, in{0} {}
};

int main() {
  const long samples{1000000l};   // total number of points in square
  trng::yarn2 r;                  // random number engine
  parallel_pi pi(r);              // functor for parallel reduce
  // parallel MC computation of pi
  tbb::parallel_reduce(tbb::blocked_range<long>(0, samples), pi, tbb::auto_partitioner());
  // print result
  std::cout << "pi = " << 4.0 * pi.in_circle() / samples << std::endl;
  return EXIT_SUCCESS;
}
