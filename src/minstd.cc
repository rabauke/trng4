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

#include <trng/minstd.hpp>

namespace trng {

  // Uniform random number generator concept
  const minstd::result_type minstd::min=1l;
  const minstd::result_type minstd::max=2147483646l;
  
  // Parameter and status classes
  
  // Equality comparable concept
  bool operator==(const minstd::status_type &S1, 
		  const minstd::status_type &S2) {
    return S1.r==S2.r;
  }

  bool operator!=(const minstd::status_type &S1, 
		  const minstd::status_type &S2) {
    return !(S1==S2);
  }
  
  // Random number engine concept
  minstd::minstd() : S() { }

  minstd::minstd(unsigned long s) : S() { 
    seed(s);
  }
  
  void minstd::seed() {
    (*this)=minstd();
  }
 
  void minstd::seed(minstd::result_type s) {
    S.r=s%2147483647ul;
    if (S.r==0)
      S.r=1;
  }
  
  // Equality comparable concept
  bool operator==(const minstd &R1, const minstd &R2) {
    return R1.S==R2.S;
  }

  bool operator!=(const minstd &R1, const minstd &R2) {
    return !(R1==R2);
  }
  
  // Other usefull methods
  const char * const minstd::name_str="minstd";
  
  const char * minstd::name() {
    return name_str;
  }
  
}

