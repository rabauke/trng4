#if __cplusplus > 201100L
// mix random number generators and distribution from
// TRNG and the C++11 standard library

#include <cstdlib>
#include <iostream>
#include <random>
#include <trng/lcg64.hpp>
#include <trng/normal_dist.hpp>

int main() {
  std::mt19937 R_cpp11;
  trng::lcg64 R_trng;
  std::normal_distribution<> N_cpp11;
  trng::normal_dist<> N_trng(0, 1);
  for (int i=0; i<10000; ++i) {
    std::cout << N_cpp11(R_cpp11) << '\t';
    std::cout << N_cpp11(R_trng) << '\t';
    std::cout << N_trng(R_cpp11) << '\t';
    std::cout << N_trng(R_trng) << '\n';
  }
  return EXIT_SUCCESS;
}
#else
#include <iostream>
int main() {
  std::cout << "Not a C++11 compiler.\n";
}
#endif
