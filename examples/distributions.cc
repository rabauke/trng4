// Copyright (c) 2000-2014, Heiko Bauke
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
#include <iterator>
#include <numeric>
#include <trng/limits.hpp>
#include <trng/math.hpp>
#include <trng/yarn2.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/uniform_int_dist.hpp>
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
#include <trng/binomial_dist.hpp>
#include <trng/negative_binomial_dist.hpp>
#include <trng/hypergeometric_dist.hpp>
#include <trng/geometric_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <trng/discrete_dist.hpp>


// integration by Simpson rule
//
template<typename iter>
typename std::iterator_traits<iter>::value_type
simpson_int(iter first, iter last) {
  typedef typename std::iterator_traits<iter>::value_type value_type;
  if (first==last)
    return value_type(0);
  value_type sum1(0), sum2(0);
  int n=std::distance(first, last);
  if (n%2==0) { // Pulcherima (3/8 rule) for even number of data points
    sum1+=(*first)*value_type(3)/value_type(8);
    ++first;
    sum1+=(*first)*value_type(9)/value_type(8);
    ++first;
    sum1+=(*first)*value_type(9)/value_type(8);
    ++first;
    sum1+=(*first)*value_type(3)/value_type(8);
    n-=3;
  }
  if (n>2) {
    sum2+=*first;
    ++first;
    --n;
    while (n>0) {
      sum2+=value_type(4)*(*first);
      ++first;
      --n;
      if (n==1) 
	sum2+=(*first);
      else
	sum2+=value_type(2)*(*first);
      ++first;
      --n;
    }
    sum2/=value_type(3);
  }
  return sum1+sum2;
}


template<typename dist>
void continous_dist_test(dist &d, const char *name) {
  typedef typename dist::result_type result_type;
  
  int samples=1024*1024*32+1;
  result_type x_min=d.icdf(0.01), x_max=d.icdf(0.99), dx=(x_max-x_min)/samples;
  std::vector<result_type> y;
  for (int i=0; i<=samples; ++i) 
    y.push_back(d.pdf(x_min+i*dx));
  result_type s=simpson_int(y.begin(), y.end())*dx;
  std::cout << "\"" << name 
	    << "\" distribution test : Int_x(0.01)^x(0.99) p(x) dx = 0.98 - "
	    << 0.98-s << '\n';
  result_type dp=result_type(1)/result_type(1024*1024);
  bool ok(true);
  for (result_type p=dp; p<1; p+=dp) {
    result_type x=d.icdf(p);
    result_type y=d.cdf(x);
    if (trng::math::abs(y-p)>16*trng::math::numeric_limits<result_type>::epsilon()) {
      ok=false;
      std::cout << x << '\t' << y << '\t' << trng::math::abs(y-p) << '\n';
      break;
    }
  }
  if (ok)
    std::cout << "\"" << name << "\" cumulative distribution test passed\n";
  else 
    std::cout << "\"" << name << "\" cumulative distribution test not passed\n"; 
}


template<typename dist>
void chi2_test(dist &d, const char *name) {
  typedef typename dist::result_type result_type;
  const int bins=128;
  const result_type dp=1./bins;
  const int N=10000;
  std::vector<result_type> quantil;
  std::vector<int> count(bins, 0);
  for (int i=1; i<bins ; ++i)
    quantil.push_back(d.icdf(dp*i));
  trng::yarn2 R;
  for (int i=0; i<N; ++i) {
    result_type x=d(R);
    int bin=0;
    while (bin<bins-1 and x>quantil[bin])
      ++bin;
    ++count[bin];
  }
  result_type c2=0;
  for (int i=0; i<bins ; ++i)
    c2+=(count[i]-N*dp)*(count[i]-N*dp)/(N*dp);
  result_type c2_p=trng::math::GammaQ(0.5*c2, 0.5*(bins-1));
  if (0.01<c2_p and c2_p<0.99)
    std::cout << "\"" << name << "\" chi-squared test passed\n";
  else {
    std::cout << "\"" << name << "\" chi-squared test not passed\n"; 
    for (int i=0; i<bins ; ++i)
      std::cout << count[i] << '\n';
  }
}


template<typename dist>
void discrete_dist_test(dist &d, const char *name) {
  bool ok=true;
  double p=0;
  int i=0;
  while (p<0.9) {
    p+=d.pdf(i);
    if (trng::math::abs(p-d.cdf(i))>16*trng::math::numeric_limits<double>::epsilon()) {
      ok=false;
      break;
    }
    ++i;
  }
  if (ok)
    std::cout << "\"" << name << "\" distribution test passed\n";
  else 
    std::cout << "\"" << name << "\" distribution test not passed\n";
}


template<typename dist>
void continous_test(dist &d, const char *name) {
  continous_dist_test(d, name);
  chi2_test(d, name);
}


template<typename dist>
void discrete_test(dist &d, const char *name) {
  discrete_dist_test(d, name);
}


int main() {
  {
    trng::uniform_dist<> d(2.0, 5.0);
    continous_test(d, "uniform distribution");
  }
  {
    trng::uniform01_dist<> d;
    continous_test(d, "uniform01 distribution");
  }
  {
    trng::exponential_dist<> d(2.0);
    continous_test(d, "exponential distribution");
  }
  {
    trng::twosided_exponential_dist<> d(2.0);
    continous_test(d, "twosided exponential distribution");
  }
  {
    trng::normal_dist<> d(5.0, 2.0);
    continous_test(d, "normal distribution");
  }
  {
    trng::truncated_normal_dist<> d(5.0, 2.0, 2.0, 6.0);
    continous_test(d, "truncated normal distribution");
  }
  {
    trng::maxwell_dist<> d(2.0);
    continous_test(d, "maxwell distribution");
  }
  {
    trng::cauchy_dist<> d(5.0, 2.0);
    continous_test(d, "cauchy distribution");
  }
  {
    trng::logistic_dist<> d(5.0, 2.0);
    continous_test(d, "logistic distribution");
  }
  {
    trng::lognormal_dist<> d(5.0, 2.0);
    continous_test(d, "lognormal distribution");
  }
  {
    trng::pareto_dist<> d(5.0, 2.0);
    continous_test(d, "pareto distribution");
  }
  {
    trng::powerlaw_dist<> d(5.0, 2.0);
    continous_test(d, "power-law distribution");
  }
  {
    trng::tent_dist<> d(5.0, 2.0);
    continous_test(d, "tent distribution");
  }
  {
    trng::extreme_value_dist<> d(5.0, 2.0);
    continous_test(d, "extreme-value distribution");
  }
  {
    trng::gamma_dist<> d(5.0, 2.0);
    continous_test(d, "gamma distribution");
  }
  {
    trng::beta_dist<> d(3.0, 2.0);
    continous_test(d, "beta distribution");
  }
  {
    trng::chi_square_dist<> d(8);
    continous_test(d, "chi-square distribution");
  }
  {
    trng::student_t_dist<> d(10);
    continous_test(d, "student-t distribution");
  }
  {
    trng::snedecor_f_dist<> d(10, 11);
    continous_test(d, "snedecor-f distribution");
  }
  {
    trng::rayleigh_dist<> d(10);
    continous_test(d, "rayleigh distribution");
  }

  {
    trng::bernoulli_dist<int> d(0.4, 0, 1);
    discrete_test(d, "bernoulli distribution");
  }
  {
    trng::binomial_dist d(0.4, 20);
    discrete_test(d, "binomial distribution");
  }
  {
    trng::negative_binomial_dist d(0.4, 20);
    discrete_test(d, "negative binomial distribution");
  }
  {
    trng::hypergeometric_dist d(10, 5, 5);
    discrete_test(d, "hypergeometric distribution");
  }
  {
    trng::geometric_dist d(0.3);
    discrete_test(d, "geometric distribution");
  }
  {
    trng::poisson_dist d(0.3);
    discrete_test(d, "poisson distribution");
  }
  return EXIT_SUCCESS;
}
