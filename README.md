# TRNG - A modern C++ pseudo random number generator library

## Key features

* fully compatible with the C++11 random number facility as defined in [`<random>`](https://en.cppreference.com/w/cpp/header/random)
* implements various pseudo random number algorithms
* supports multiple streams of random numbers for parallel (multi-threaded) applications
* does not depend on a specific parallelization technique, may be used with any threading library or MPI
* pseudo random numbers can be sampled from different distributions
* bindings for the R programming language provided via [rTRNG package](https://cran.r-project.org/web/packages/rTRNG/index.html), see also [this blog post](https://mirai-solutions.ch/news/2019/06/10/rTRNG-avanced-parallel-RNG-R/) or [this presentation](https://docs.microsoft.com/en-US/events/user-international-r-user-conferences-user-international-r-user-2017-conference/rtrng-advanced-parallel-random-number-generation-in-r)


## Documentation

For installation instructions and documentation read trng.pdf in the
doc directory. 
