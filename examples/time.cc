// Copyright (C) 2006 Heiko Bauke <heiko.bauke@physik.uni-magdeburg.de>
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

#include <cstdlib>
#include <iostream>
#include <exception>
#include <string>
#include <trng/config.hpp>
#include <trng/lcg64.hpp>
#include <trng/mrg2.hpp>
#include <trng/mrg3.hpp>
#include <trng/mrg3s.hpp>
#include <trng/mrg4.hpp>
#include <trng/mrg5.hpp>
#include <trng/mrg5s.hpp>
#include <trng/yarn2.hpp>
#include <trng/yarn3.hpp>
#include <trng/yarn3s.hpp>
#include <trng/yarn4.hpp>
#include <trng/yarn5.hpp>
#include <trng/yarn5s.hpp>
#include <trng/lagfib4xor.hpp>

#if defined HAVE_BOOST
# include <boost/random/linear_congruential.hpp>
# include <boost/random/additive_combine.hpp>
# include <boost/random/inversive_congruential.hpp>
# include <boost/random/mersenne_twister.hpp>
# include <boost/random/lagged_fibonacci.hpp>
# include <boost/random/shuffle_output.hpp>
#endif

#if defined __unix__
# include <unistd.h>
# include <sys/time.h>
# include <sys/times.h>
#else
# include <ctime>
#endif

class timer {
private:
  const double _resolution;
  double _t;
  double get_time() {
#if defined __unix__
    struct timeval  tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return static_cast<double>(tv.tv_sec)+static_cast<double>(tv.tv_usec)*1e-6;
#else
    return static_cast<double>(std::clock())*_resolution;
#endif
  }
public:
  void reset() { _t=get_time(); }
  double time() { return get_time()-_t; }
  double resolution() const { return _resolution; }
  timer() :
#if defined __unix__
    _resolution(1e-6),
#else
    _resolution(1.0/CLOCKS_PER_SEC),
#endif
    _t(get_time()) { }
};

template<typename R>
typename R::result_type time_main(R &r, std::string name) {
  typename R::result_type s(0);
  timer T;
  long max((1<<26));
  for (long i(0); i<max; ++i)
    s+=r();
  while (name.length()<20)
    name+=' ';
  std::cout << name << '\t'
	    << 1e-6*max/T.time() << '\n';
  return s;
}

int main() {
  std::cout << "Generator               10^6 random numbers per second\n"
	    << "======================================================\n";
  try {
    { trng::lcg64    r;  time_main(r, "trng::lcg64"); }
    { trng::mrg2     r;  time_main(r, "trng::mrg2"); }
    { trng::mrg3     r;  time_main(r, "trng::mrg3"); }
    { trng::mrg3s    r;  time_main(r, "trng::mrg3s"); }
    { trng::mrg4     r;  time_main(r, "trng::mrg4"); }
    { trng::mrg5     r;  time_main(r, "trng::mrg5"); }
    { trng::mrg5s    r;  time_main(r, "trng::mrg5s"); }
    { trng::yarn2    r;  time_main(r, "trng::yarn2"); }
    { trng::yarn3    r;  time_main(r, "trng::yarn3"); }
    { trng::yarn3s   r;  time_main(r, "trng::yarn3s"); }
    { trng::yarn4    r;  time_main(r, "trng::yarn4"); }
    { trng::yarn5    r;  time_main(r, "trng::yarn5"); }
    { trng::yarn5s   r;  time_main(r, "trng::yarn5s"); }
    { trng::Ziff_ul  r;  time_main(r, "trng::Ziff_ul"); }
#if defined HAVE_BOOST
    { boost::minstd_rand    r; time_main(r, "boost::minstd_rand"); }
    { boost::ecuyer1988     r; time_main(r, "boost::ecuyer1988"); }
    { boost::kreutzer1986   r; time_main(r, "boost::kreutzer1986"); }
    { boost::hellekalek1995 r; time_main(r, "boost::hellekalek1995"); }
    { boost::mt11213b       r; time_main(r, "boost::mt11213b"); }
    { boost::mt19937        r; time_main(r, "boost::mt19937"); }
    { boost::lagged_fibonacci607   r; time_main(r, "boost::lagged_fibonacci607"); }
    { boost::lagged_fibonacci1279  r; time_main(r, "boost::lagged_fibonacci1279"); }
    { boost::lagged_fibonacci2281  r; time_main(r, "boost::lagged_fibonacci2281"); }
    { boost::lagged_fibonacci3217  r; time_main(r, "boost::lagged_fibonacci3217"); }
    { boost::lagged_fibonacci4423  r; time_main(r, "boost::lagged_fibonacci4423"); }
    { boost::lagged_fibonacci9689  r; time_main(r, "boost::lagged_fibonacci9689"); }
    { boost::lagged_fibonacci19937 r; time_main(r, "boost::lagged_fibonacci19937"); }
    { boost::lagged_fibonacci23209 r; time_main(r, "boost::lagged_fibonacci23209"); }
    { boost::lagged_fibonacci44497 r; time_main(r, "boost::lagged_fibonacci44497"); }
#endif
  } 
  catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
  }
  catch (...) {
    std::cerr << "something went wrong" << std::endl;
  }
  return EXIT_SUCCESS;
}


