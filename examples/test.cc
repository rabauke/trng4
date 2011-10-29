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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
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

template<typename R> 
void test() {
  R engine;
  std::string fname(R::name());
  fname+=".dat";
  std::ofstream out(fname.c_str());
  for (int i=0; i<1000; ++i)
    out << engine() << '\n';
}

int main() {
  test<trng::lcg64>();
  test<trng::lcg64_shift>();
  test<trng::mrg2>();
  test<trng::mrg3>();
  test<trng::mrg3s>();
  test<trng::mrg4>();
  test<trng::mrg5>();
  test<trng::mrg5s>();
  test<trng::yarn2>();
  test<trng::yarn3>();
  test<trng::yarn3s>();
  test<trng::yarn4>();
  test<trng::yarn5>();
  test<trng::yarn5s>();
  return EXIT_SUCCESS;
}
