#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <trng/lcg64.hpp>
#include <trng/discrete_dist.hpp>

int main() {
  std::vector<double> p;  // stores relative probabilities
  // populate vector with relative probabilities
  p.push_back(1);
  p.push_back(3.25);
  p.push_back(5);
  p.push_back(6.5);
  p.push_back(7);
  p.push_back(2);
  // discrete distribution object
  trng::discrete_dist dist(p.begin(), p.end());
  // random number generator
  trng::lcg64 r;
  // draw and some random numbers 
  std::vector<int> count(p.size(), 0);
  const int samples=10000;
  for (int i=0; i<samples; ++i) {
    int x=dist(r);  // draw a random number
    ++count[x];     // count
  }
  // print results
  std::cout << "value\t\tprobability\tcount\t\tempirical probability\n"
	    << "=====\t\t===========\t=====\t\t=====================\n";
  for (std::vector<int>::size_type i=0; i<count.size(); ++i) {
    std::cout << std::setprecision(3) 
	      << i << "\t\t"
	      << dist.pdf(i) << "\t\t"
	      << count[i] << "\t\t"
	      << static_cast<double>(count[i])/samples << '\n';
  }
  return EXIT_SUCCESS;
}
