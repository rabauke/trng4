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
#include <iomanip>
#include <trng/lcg64.hpp>
#include <trng/discrete_dist.hpp>

int main() {
  const int n{6};
  // stores relative probabilities
  const double p[n]{1., 3.25, 5., 6.5, 7., 2.};
  // discrete distribution object
  trng::discrete_dist dist(std::begin(p), std::end(p));
  // random number generator
  trng::lcg64 r;
  // draw some random numbers
  int count[n]{0, 0, 0, 0, 0, 0};
  const int samples{10000};
  for (int i{0}; i < samples; ++i) {
    const int x{dist(r)};  // draw a random number
    ++count[x];            // count
  }
  // print results
  std::cout << "value\t\tprobability\tcount\t\tempirical probability\n"
            << "=====\t\t===========\t=====\t\t=====================\n";
  for (int i{0}; i < n; ++i) {
    std::cout << std::setprecision(3) << i << "\t\t" << dist.pdf(i) << "\t\t" << count[i]
              << "\t\t" << static_cast<double>(count[i]) / samples << '\n';
  }
  return EXIT_SUCCESS;
}
