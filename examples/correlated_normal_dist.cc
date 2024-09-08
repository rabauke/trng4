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
#include <iomanip>
#include <vector>
#include <trng/lcg64.hpp>
#include <trng/correlated_normal_dist.hpp>

double covariance(const std::vector<double>& v1, const std::vector<double>& v2);

double covariance(const std::vector<double>& v1, const std::vector<double>& v2) {
  const std::vector<double>::size_type n{v1.size()};
  double m1{0.0}, m2{0.0}, c{0.0};
  for (std::vector<double>::size_type i{0}; i < n; ++i) {
    m1 += v1[i] / n;
    m2 += v2[i] / n;
  }
  for (std::vector<double>::size_type i{0}; i < n; ++i)
    c += (v1[i] - m1) * (v2[i] - m2) / n;
  return c;
}

int main() {
  const int d{4};
  // covariance matrix
  const double sig[d][d]{{2.0, -0.5, 0.3, -0.3},
                         {-0.5, 3.0, -0.3, 0.3},
                         {0.3, -0.3, 1.0, -0.3},
                         {-0.3, 0.3, -0.3, 1.0}};
  trng::correlated_normal_dist<> D(&sig[0][0], &sig[d - 1][d - 1] + 1);
  trng::lcg64 R;

  std::vector<double> x1, x2, x3, x4;
  // generate 4-tuples of correlated normal variables
  for (int i{0}; i < 1000000; ++i) {
    x1.push_back(D(R));
    x2.push_back(D(R));
    x3.push_back(D(R));
    x4.push_back(D(R));
  }
  // print (empirical) covariance matrix
  std::cout << std::setprecision(4) << covariance(x1, x1) << '\t' << covariance(x1, x2) << '\t'
            << covariance(x1, x3) << '\t' << covariance(x1, x4) << '\n'
            << covariance(x2, x1) << '\t' << covariance(x2, x2) << '\t' << covariance(x2, x3)
            << '\t' << covariance(x2, x4) << '\n'
            << covariance(x3, x1) << '\t' << covariance(x3, x2) << '\t' << covariance(x3, x3)
            << '\t' << covariance(x3, x4) << '\n'
            << covariance(x4, x1) << '\t' << covariance(x4, x2) << '\t' << covariance(x4, x3)
            << '\t' << covariance(x4, x4) << '\n';
  return EXIT_SUCCESS;
}
