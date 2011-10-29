// Copyright (C) 2001-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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
#include <vector>
#include <algorithm>
#include <trng/config.hpp>
#include <trng/yarn2.hpp>
#include <trng/uniform_int_dist.hpp>
#if defined TRNG_HAVE_BOOST
  #include <boost/bind.hpp>
#else

  // helper class 
  template<typename PRN_dist_t, typename PRN_engine_t>
  class binder_cl {
    PRN_dist_t &dist;
    PRN_engine_t &engine;
  public:
    binder_cl(PRN_dist_t &dist, PRN_engine_t &engine) : dist(dist), engine(engine) {
    }
    typename PRN_dist_t::result_type operator()() {
      return dist(engine);
    }
  };

  // convenience function
  template<typename PRN_dist_t, typename PRN_engine_t>
  inline 
  binder_cl<PRN_dist_t, PRN_engine_t> make_binder(PRN_dist_t &dist, PRN_engine_t &engine) {
    return binder_cl<PRN_dist_t, PRN_engine_t>(dist, engine);
  }

#endif


// print an iterator range to stdout
template<typename iter>
void print_range(iter i1, iter i2) {
  while (i1!=i2) std::cout << (*(i1++)) << '\t';
  std::cout << "\n\n";
}

int main() {
  trng::yarn2 R;
  trng::uniform_int_dist U(0, 100);
  std::vector<long> v(10);
  
  std::cout << "random number generation by call operator\n";
  for (std::vector<long>::size_type i=0; i<v.size(); ++i)
    v[i]=U(R);
  print_range(v.begin(), v.end());
  std::vector<long> w(12);
#if defined TRNG_HAVE_BOOST
  std::cout << "random number generation by std::generate\n";
  std::generate(w.begin(), w.end(), boost::bind(U, boost::ref(R)));
  print_range(w.begin(), w.end());
  std::cout << "random number generation by std::generate\n";
  std::generate(w.begin(), w.end(), boost::bind(U, boost::ref(R)));
  print_range(w.begin(), w.end());
#else
  std::cout << "random number generation by std::generate\n";
  std::generate(w.begin(), w.end(), make_binder(U, R));
  print_range(w.begin(), w.end());
  std::cout << "random number generation by std::generate\n";
  std::generate(w.begin(), w.end(), make_binder(U, R));
  print_range(w.begin(), w.end());
#endif
  std::cout << "same sequence as above, but in a random shuffled order\n";
  std::random_shuffle(w.begin(), w.end(), R);
  print_range(w.begin(), w.end());
  return EXIT_SUCCESS;
}
