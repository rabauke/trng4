# TRNG - A modern C++ pseudo random number generator library

## Key features

* fully compatible with the C++11 random number facility as defined in [`<random>`](https://en.cppreference.com/w/cpp/header/random)
* implements various pseudo random number algorithms
* supports multiple streams of random numbers for parallel (multithreaded) applications
* does not depend on a specific parallelization technique, may be used with any threading library or MPI
* pseudo random numbers can be sampled from many different discrete and continuous distributions
* employs the CMake build system and features [CMake package](https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html) support
* bindings for the R programming language provided via [rTRNG package](https://cran.r-project.org/web/packages/rTRNG/index.html), see also [this blog post](https://mirai-solutions.ch/news/2019/06/10/rTRNG-avanced-parallel-RNG-R/) or [this presentation](https://docs.microsoft.com/en-US/events/user-international-r-user-conferences-user-international-r-user-2017-conference/rtrng-advanced-parallel-random-number-generation-in-r)

## Example

TRNG classes can be used as a drop-in replacement for classes declared in the `random` header 
 file of the C++ standard library.  In addition, the TRNG random number generators provide `jump` and `split` methods for constructing independent streams of pseudo random numbers for parallel Monte Carlo simulations.  The following code illustrates the use of TRNG pseudo random number generators for the parallel [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method) calculation of pi.

```c++
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

int main() {
  const long samples{1000000l};  // total number of points in square
  long in{0l};                   // number of points in circle
  // distribute workload over all processes and make a global reduction
#pragma omp parallel reduction(+ : in) default(none)
  {
    trng::yarn2 r;                          // random number engine
    const int size{omp_get_num_threads()};  // get total number of processes
    const int rank{omp_get_thread_num()};   // get rank of current process
    trng::uniform01_dist<> u;               // random number distribution
    r.jump(2 * (rank * samples / size));    // jump ahead
    // throw random points into square
    for (long i{rank * samples / size}; i < (rank + 1) * samples / size; ++i) {
      const double x{u(r)}, y{u(r)};  // choose random x- and y-coordinates
      if (x * x + y * y <= 1.0)       // is point in circle?
        ++in;                         // increase thread-local counter
    }
  }
  // print result
  std::cout << "pi = " << 4.0 * in / samples << std::endl;
  return EXIT_SUCCESS;
}
```

## Documentation

For installation instructions and further documentation see the trng.pdf file in the
doc directory. 
