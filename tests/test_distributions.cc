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

#include <vector>
#include <iterator>
#include <limits>
#include <cmath>
#include <numeric>
#include <ciso646>
#include <sstream>

#include <catch2/catch.hpp>

#include <trng/uniform_dist.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/special_functions.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/exponential_dist.hpp>
#include <trng/twosided_exponential_dist.hpp>
#include <trng/normal_dist.hpp>
#include <trng/truncated_normal_dist.hpp>
#include <trng/maxwell_dist.hpp>
#include <trng/cauchy_dist.hpp>
#include <trng/logistic_dist.hpp>
#include <trng/lognormal_dist.hpp>
#include <trng/pareto_dist.hpp>
#include <trng/powerlaw_dist.hpp>
#include <trng/tent_dist.hpp>
#include <trng/weibull_dist.hpp>
#include <trng/extreme_value_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <trng/beta_dist.hpp>
#include <trng/chi_square_dist.hpp>
#include <trng/student_t_dist.hpp>
#include <trng/snedecor_f_dist.hpp>
#include <trng/rayleigh_dist.hpp>
#include <trng/bernoulli_dist.hpp>
#include <trng/uniform_int_dist.hpp>
#include <trng/binomial_dist.hpp>
#include <trng/negative_binomial_dist.hpp>
#include <trng/hypergeometric_dist.hpp>
#include <trng/geometric_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <trng/zero_truncated_poisson_dist.hpp>
#include <trng/discrete_dist.hpp>


// integration by Simpson rule
//
template<typename iter>
typename std::iterator_traits<iter>::value_type simpson_int(iter first, iter last) {
  typedef typename std::iterator_traits<iter>::value_type value_type;
  if (first == last)
    return value_type{0};
  value_type sum{0};
  auto n{std::distance(first, last)};
  if (n % 2 == 0) {  // Pulcherima (3/8 rule) for even number of data points
    sum += (*first) * value_type(3) / value_type(8);
    ++first;
    sum += (*first) * value_type(9) / value_type(8);
    ++first;
    sum += (*first) * value_type(9) / value_type(8);
    ++first;
    sum += (*first) * value_type(3) / value_type(8);
    n -= 3;
  }
  if (n > 2) {
    sum *= value_type(3);
    sum += *first;
    ++first;
    --n;
    while (n > 0) {
      sum += value_type(4) * (*first);
      ++first;
      --n;
      if (n == 1)
        sum += (*first);
      else
        sum += value_type(2) * (*first);
      ++first;
      --n;
    }
    sum /= value_type(3);
  }
  return sum;
}


// p value of Chi-squared test
double chi_percentil(const std::vector<double> &p, const std::vector<int> &count) {
  using size_type = std::vector<double>::size_type;
  const int N{std::accumulate(count.begin(), count.end(), 0)};
  const size_type bins{p.size()};
  double c2{0};
  for (size_type i{0}; i < bins; ++i) {
    const double expected{N * p[i]};
    c2 += (count[i] - expected) * (count[i] - expected) / expected;
  }
  const double c2_p{trng::math::GammaQ(c2 / 2, static_cast<double>(bins - 1) / 2)};
  return c2_p;
}


template<typename dist>
void continuous_dist_test_integrate_pdf(const dist &d) {
  using result_type = typename dist::result_type;
  // const int samples{1024 * 16 + 1};
  const int samples(static_cast<int>(std::min(
      1024ll * 1024, static_cast<long long>(std::round(
                         1 / std::sqrt(std::numeric_limits<result_type>::epsilon()))))));
  const result_type x_min{d.icdf(result_type(1) / result_type(100))};
  const result_type x_max{d.icdf(result_type(99) / result_type(100))};
  const result_type dx{(x_max - x_min) / samples};
  std::vector<result_type> y;
  for (int i{0}; i <= samples; ++i)
    y.push_back(d.pdf(x_min + i * dx));
  const result_type s{simpson_int(y.begin(), y.end()) * dx};
  const result_type tol{result_type(64) / result_type(samples) / result_type(samples)};
  REQUIRE(std::abs(s - result_type(98) / result_type(100)) < tol);
}


template<typename dist>
void continuous_dist_test_icdf(const dist &d) {
  using result_type = typename dist::result_type;
  const int bins{1024 * 1024};
  const result_type dp{result_type(1) / result_type(bins)};
  const auto eps{256 * std::numeric_limits<result_type>::epsilon()};
  for (int i{1}; i < bins; ++i) {
    const result_type p{i * dp};
    const result_type x{d.icdf(p)};
    const result_type y{d.cdf(x)};
    if (not(std::abs(y - p) < eps) or i == 1)
      REQUIRE(std::abs(y - p) < eps);
  }
}


// failure of the chi2 test does not necessarily imply an error, may happen just by chance
template<typename dist>
void continuous_dist_test_chi2_test(dist &d) {
  using result_type = typename dist::result_type;
  const int bins{128};
  const result_type dp{result_type(1) / result_type(bins)};
  const int N{10000};
  std::vector<result_type> quantil;
  std::vector<int> count(bins, 0);
  for (int i{1}; i < bins; ++i)
    quantil.push_back(d.icdf(dp * i));
  trng::lcg64_shift R;
  for (int i{0}; i < N; ++i) {
    const result_type x{d(R)};
    int bin{0};
    while (bin < bins - 1 and x > quantil[bin])
      ++bin;
    ++count[bin];
  }
  const std::vector<double> p(bins, 1.0 / bins);
  const double c2_p{chi_percentil(p, count)};
  REQUIRE((0.01 < c2_p and c2_p < 0.99));
}


template<typename dist>
void continuous_dist_test_streamable(dist &d) {
  typename dist::param_type p_new;
  std::stringstream str;
  str << d.param();
  str >> p_new;
  dist d_new{p_new};
  REQUIRE((d == d_new));
}


template<typename T>
void continuous_dist_test(T &d) {
  SECTION("integrate pdf") { continuous_dist_test_integrate_pdf(d); }
  SECTION("icdf") { continuous_dist_test_icdf(d); }
  SECTION("chi2 test") { continuous_dist_test_chi2_test(d); }
  SECTION("streamable") { continuous_dist_test_streamable(d); }
}


template<typename dist>
void discrete_dist_test_pdf(dist &d) {
  int i{d.min()};
  double P{d.cdf(i)};
  while (P < 0.95) {
    double p{P};
    if (i > d.min())
      p -= d.cdf(i - 1);
    const double diff{p - d.pdf(i)};
    REQUIRE(diff < 128 * std::numeric_limits<double>::epsilon());
    ++i;
    P = d.cdf(i);
  }
}


// failure of the chi2 test does not necessarily imply an error, may happen just by chance
template<typename dist>
void discrete_dist_test_chi2_test(dist &d) {
  using result_type = typename dist::result_type;
  std::vector<double> p;
  {
    int i{d.min()};
    while (i <= d.max()) {
      p.push_back(d.pdf(i));
      double P{d.cdf(i)};
      if (P > 0.99)
        break;
      ++i;
    }
    if (i < d.max())
      p.push_back(1.0 - d.cdf(i));
  }
  const int bins{static_cast<int>(p.size())};
  const int N{10000};
  std::vector<int> count(bins, 0);
  trng::lcg64_shift R(100ull);
  for (int i{0}; i < N; ++i) {
    const result_type x{d(R) - d.min()};
    int bin{std::min(x, bins - 1)};
    ++count[bin];
  }
  // merge bins with very small count numbers
  while (p.size() > 2) {
    auto min{std::min_element(count.begin(), count.end())};
    if (*min > 8)
      break;
    auto pos{min - count.begin()};
    const double p_old{p[pos]};
    const int count_old{count[pos]};
    p.erase(p.begin() + pos);
    count.erase(count.begin() + pos);
    min = std::min_element(count.begin(), count.end());
    pos = min - count.begin();
    p[pos] += p_old;
    count[pos] += count_old;
  }
  const double c2_p{chi_percentil(p, count)};
  REQUIRE((0.01 < c2_p and c2_p < 0.99));
}


template<typename dist>
void discrete_dist_test_streamable(dist &d) {
  typename dist::param_type p_new;
  std::stringstream str;
  str << d.param();
  str >> p_new;
  dist d_new{p_new};
  REQUIRE(d == d_new);
}


template<typename T>
void discrete_dist_test(T &d) {
  SECTION("test pdf & cdf") { discrete_dist_test_pdf(d); }
  SECTION("chi2_test") { discrete_dist_test_chi2_test(d); }
  SECTION("streamable") { discrete_dist_test_streamable(d); }
}


TEMPLATE_TEST_CASE("continuous distributions", "", float, double, long double) {
  SECTION("uniform_dist") {
    trng::uniform_dist<TestType> d(TestType(2), TestType(5));
    continuous_dist_test(d);
  }

  SECTION("uniform01_dist") {
    trng::uniform01_dist<TestType> d;
    continuous_dist_test(d);
  }

  SECTION("exponential_dist") {
    trng::exponential_dist<TestType> d(TestType(2));
    continuous_dist_test(d);
  }

  SECTION("twosided_exponential_dist") {
    trng::twosided_exponential_dist<TestType> d(TestType(2));
    continuous_dist_test(d);
  }

  SECTION("normal_dist") {
    trng::normal_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("truncated_normal_dist") {
    trng::truncated_normal_dist<TestType> d(TestType(5), TestType(2), TestType(2), TestType(6));
    continuous_dist_test(d);
  }

  SECTION("maxwell_dist") {
    trng::maxwell_dist<TestType> d(TestType(2));
    continuous_dist_test(d);
  }

  SECTION("cauchy_dist") {
    trng::cauchy_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("logistic_dist") {
    trng::logistic_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("lognormal_dist") {
    trng::lognormal_dist<TestType> d(TestType(1), TestType(1) / TestType(2));
    continuous_dist_test(d);
  }

  SECTION("pareto_dist") {
    trng::pareto_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("powerlaw_dist") {
    trng::powerlaw_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("tent_dist") {
    trng::tent_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("weibull_dist") {
    trng::weibull_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("extreme_value_dist") {
    trng::extreme_value_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("gamma_dist") {
    trng::gamma_dist<TestType> d(TestType(5), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("beta_dist") {
    trng::beta_dist<TestType> d(TestType(3), TestType(2));
    continuous_dist_test(d);
  }

  SECTION("chi_square_dist") {
    trng::chi_square_dist<TestType> d(38);
    continuous_dist_test(d);
  }

  SECTION("student_t_dist") {
    trng::student_t_dist<TestType> d(10);
    continuous_dist_test(d);
  }

  SECTION("snedecor_f_dist") {
    trng::snedecor_f_dist<TestType> d(10, 11);
    continuous_dist_test(d);
  }

  SECTION("rayleigh_dist") {
    trng::rayleigh_dist<TestType> d(TestType(10));
    continuous_dist_test(d);
  }
}


TEMPLATE_TEST_CASE("discrete distributions", "", float, double, long double) {
  SECTION("bernoulli_dist") {
    trng::bernoulli_dist<int> d(0.4);
    discrete_dist_test(d);
  }

  SECTION("uniform_int_dist") {
    trng::uniform_int_dist d(8, 100);
    discrete_dist_test(d);
  }

  SECTION("binomial_dist") {
    trng::binomial_dist d(0.4, 20);
    discrete_dist_test(d);
  }

  SECTION("negative_binomial_dist") {
    trng::negative_binomial_dist d(0.4, 20);
    discrete_dist_test(d);
  }

  SECTION("hypergeometric_dist") {
    trng::hypergeometric_dist d(10, 5, 5);
    discrete_dist_test(d);
  }

  SECTION("geometric_dist") {
    trng::geometric_dist d(0.3);
    discrete_dist_test(d);
  }

  SECTION("poisson_dist") {
    trng::poisson_dist d(2.125);
    discrete_dist_test(d);
  }

  SECTION("zero_truncated_poisson_dist") {
    trng::zero_truncated_poisson_dist d(2.125);
    discrete_dist_test(d);
  }

  SECTION("discrete_dist") {
    std::vector<int> p{1, 2, 3, 4, 5, 4, 3, 2, 1};
    trng::discrete_dist d(p.begin(), p.end());
    discrete_dist_test(d);
  }
}
