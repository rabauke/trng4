// Copyright (c) 2000-2021, Heiko Bauke
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
#include <vector>
#include <trng/lcg64.hpp>
#include <trng/bernoulli_dist.hpp>

typedef enum { head = 0, tail = 1 } coin;

int main() {
  // discrete distribution object
  trng::bernoulli_dist<coin> biased_coin(0.51, head, tail);
  // random number generator
  trng::lcg64 r;
  // draw some random numbers
  std::vector<int> count(2, 0);
  const int samples{100000};
  for (int i = 0; i < samples; ++i) {
    const int x{biased_coin(r)};  // draw a random number
    ++count[x];                   // count
  }
  // print results
  std::cout << "value\t\tprobability\tcount\t\tempirical probability\n"
            << "=====\t\t===========\t=====\t\t=====================\n";
  for (std::vector<int>::size_type i = 0; i < count.size(); ++i)
    std::cout << std::setprecision(3) << i << "\t\t" << biased_coin.pdf(static_cast<coin>(i))
              << "\t\t" << count[i] << "\t\t" << static_cast<double>(count[i]) / samples
              << '\n';
  return EXIT_SUCCESS;
}
