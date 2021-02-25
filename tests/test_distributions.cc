// Copyright (c) 2000-2020, Heiko Bauke
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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

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
#include <trng/discrete_dist.hpp>


using floats = boost::mpl::list<float, double, long double>;

// integration by Simpson rule
//
template<typename iter>
typename std::iterator_traits<iter>::value_type simpson_int(iter first, iter last) {
  typedef typename std::iterator_traits<iter>::value_type value_type;
  if (first == last)
    return value_type{0};
  value_type sum{0};
  auto n = std::distance(first, last);
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
bool continous_dist_test_integrate_pdf(const dist &d) {
  using result_type = typename dist::result_type;
  // const int samples{1024 * 16 + 1};
  const int samples(static_cast<int>(std::min(1024ll * 1024,
                             static_cast<long long>(std::round(
                         1 / std::sqrt(std::numeric_limits<result_type>::epsilon()))))));
  const result_type x_min{d.icdf(result_type(1) / result_type(100))};
  const result_type x_max{d.icdf(result_type(99) / result_type(100))};
  const result_type dx{(x_max - x_min) / samples};
  std::vector<result_type> y;
  for (int i{0}; i <= samples; ++i)
    y.push_back(d.pdf(x_min + i * dx));
  const result_type s{simpson_int(y.begin(), y.end()) * dx};
  const result_type tol{result_type(64) / result_type(samples) / result_type(samples)};
  return std::abs(s - result_type(98) / result_type(100)) < tol;
}


template<typename dist>
boost::test_tools::predicate_result continous_dist_test_icdf(const dist &d) {
  using result_type = typename dist::result_type;
  const int bins{1024 * 1024};
  const result_type dp{result_type(1) / result_type(bins)};
  for (int i{1}; i < bins; ++i) {
    const result_type p{i * dp};
    const result_type x{d.icdf(p)};
    const result_type y{d.cdf(x)};
    if (std::abs(y - p) > 256 * std::numeric_limits<result_type>::epsilon()) {
      boost::test_tools::predicate_result res(false);
      res.message() << "cdf(icdf(p)) != p  for  p = " << p << " with cdf(icdf(p)) = " << y
                    << " and |cdf(icdf(p)) - p| = " << std::abs(y - p);
      return res;
    }
  }
  return true;
}


// failure of the chi2 test does not nessiarily imply an error, may happen just by chance
template<typename dist>
bool continous_dist_test_chi2_test(dist &d) {
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
  return 0.01 < c2_p and c2_p < 0.99;
}


template<typename dist>
boost::test_tools::predicate_result discrete_dist_test(dist &d) {
  int i{d.min()};
  double P{d.cdf(i)};
  while (P < 0.95) {
    double p{P};
    if (i > d.min())
      p -= d.cdf(i - 1);
    const double diff{p - d.pdf(i)};
    if (diff > 128 * std::numeric_limits<double>::epsilon()) {
      boost::test_tools::predicate_result res(false);
      res.message() << "insufficient accuracy, i = " << i << ", pdf(i) = " << d.pdf(i)
                    << ", cdf(i) - cdf(i-1) = " << p
                    << ", |cdf(i) - cdf(i-1) - pdf(i)| = " << diff;
      return res;
    }
    ++i;
    P = d.cdf(i);
  }
  return true;
}


// failure of the chi2 test does not nessiarily imply an error, may happen just by chance
template<typename dist>
bool discrete_dist_test_chi2_test(dist &d) {
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
    auto min = std::min_element(count.begin(), count.end());
    if (*min > 8)
      break;
    auto pos = min - count.begin();
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
  return 0.01 < c2_p and c2_p < 0.99;
}

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_distributions)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_continous_distributions)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_uniform_dist, T, floats) {
  trng::uniform_dist<T> d(T(2), T(5));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_uniform01_dist, T, floats) {
  trng::uniform01_dist<T> d;
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_exponential_dist, T, floats) {
  trng::exponential_dist<T> d(T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_twosided_exponential_dist, T, floats) {
  trng::twosided_exponential_dist<T> d(T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_normal_dist, T, floats) {
  trng::normal_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_truncated_normal_dist, T, floats) {
  trng::truncated_normal_dist<T> d(T(5), T(2), T(2), T(6));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_maxwell_dist, T, floats) {
  trng::maxwell_dist<T> d(T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_cauchy_dist, T, floats) {
  trng::cauchy_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_logistic_dist, T, floats) {
  trng::logistic_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_lognormal_dist, T, floats) {
  trng::lognormal_dist<T> d(T(1), T(1) / T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_pareto_dist, T, floats) {
  trng::pareto_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_powerlaw_dist, T, floats) {
  trng::powerlaw_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_tent_dist, T, floats) {
  trng::tent_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_weibull_dist, T, floats) {
  trng::weibull_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_extreme_value_dist, T, floats) {
  trng::extreme_value_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_gamma_dist, T, floats) {
  trng::gamma_dist<T> d(T(5), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta_dist, T, floats) {
  trng::beta_dist<T> d(T(3), T(2));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_chi_square_dist, T, floats) {
  trng::chi_square_dist<T> d(38);
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_student_t_dist, T, floats) {
  trng::student_t_dist<T> d(10);
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_snedecor_f_dist, T, floats) {
  trng::snedecor_f_dist<T> d(10, 11);
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_rayleigh_dist, T, floats) {
  trng::rayleigh_dist<T> d(T(10));
  BOOST_TEST(continous_dist_test_integrate_pdf(d));
  BOOST_TEST(continous_dist_test_icdf(d));
  BOOST_TEST(continous_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_discrete_distributions)

BOOST_AUTO_TEST_CASE(test_bernoulli_dist) {
  trng::bernoulli_dist<int> d(0.4);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_uniform_int_dist) {
  trng::uniform_int_dist d(8, 100);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_binomial_dist) {
  trng::binomial_dist d(0.4, 20);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_negative_binomial_dist) {
  trng::negative_binomial_dist d(0.4, 20);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_hypergeometric_dist) {
  trng::hypergeometric_dist d(10, 5, 5);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_geometric_dist) {
  trng::geometric_dist d(0.3);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_poisson_dist) {
  trng::poisson_dist d(0.3);
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_CASE(test_discrete_dist) {
  std::vector<int> p{1, 2, 3, 4, 5, 4, 3, 2, 1};
  trng::discrete_dist d(p.begin(), p.end());
  BOOST_TEST(discrete_dist_test(d));
  BOOST_TEST(discrete_dist_test_chi2_test(d));
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
