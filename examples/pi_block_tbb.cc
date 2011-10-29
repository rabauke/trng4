// Copyright (C) 2001-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License in
// version 2 as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
//

#include <trng/config.hpp>
#if defined TRNG_HAVE_TBB

#include <cstdlib>
#include <iostream>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

class parallel_pi {
  trng::uniform01_dist<> u;  // random number distribution
  long in;
  const trng::yarn2 &r;
public:
  void operator()(const tbb::blocked_range<long> &range) {
    trng::yarn2 r_local(r);                // local copy of random number engine
    r_local.jump(2*range.begin());         // jump ahead
    for (long i=range.begin(); i!=range.end(); ++i) {
      double x=u(r_local), y=u(r_local);   // choose random x- and y-coordinates
      if (x*x+y*y<=1.0)                    // is point in circle?
	++in;                              // increase thread-local counter
    }
  }
  // join threds and counters
  void join(const parallel_pi &other) {  
    in+=other.in;  
  }
  long in_circle() const {  
    return in;  
  }
  parallel_pi(const trng::yarn2 &r) : r(r), in(0) {  
  }
  parallel_pi(const parallel_pi &other, tbb::split) : r(other.r), in(0) {  
  }
};

int main(int argc, char *argv[]) {
  tbb::task_scheduler_init init;        // initiallize TBB task scheduler
  const long samples=1000000l;          // total number of points in square
  trng::yarn2 r;                        // random number engine
  parallel_pi pi(r);                    // functor for parallel reduce
  // parallel MC computation of pi 
  tbb::parallel_reduce(tbb::blocked_range<long>(0, samples), pi, tbb::auto_partitioner());
  // print result
  std::cout << "pi = " << 4.0*pi.in_circle()/samples << std::endl;
  return EXIT_SUCCESS;
}

#else

#include <cstdlib>
#include <iostream>

int main() {
  std::cerr << "Sorry, Intel Threading Buildings Blocks are not installed on your system.\n";
  return EXIT_FAILURE;
}

#endif
