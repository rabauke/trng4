// Copyright (C) 2000-2008 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#include <trng/yarn3.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // EqualityComparable concept
  bool operator==(const yarn3::parameter_type &P1,
                  const yarn3::parameter_type &P2) {
    return P1.a1==P2.a1 and P1.a2==P2.a2 and P1.a3==P2.a3;
  }

  bool operator!=(const yarn3::parameter_type &P1,
                  const yarn3::parameter_type &P2) {
    return not (P1==P2);
  }

  // Equality comparable concept
  bool operator==(const yarn3::status_type &S1,
                  const yarn3::status_type &S2) {
    return S1.r1==S2.r1 and S1.r2==S2.r2 and S1.r3==S2.r3;
  }

  bool operator!=(const yarn3::status_type &S1,
                  const yarn3::status_type &S2) {
    return not (S1==S2);
  }

  const yarn3::parameter_type
  yarn3::LEcuyer1=parameter_type(2021422057l, 1826992351l, 1977753457l);
  const yarn3::parameter_type
  yarn3::LEcuyer2=parameter_type(1476728729l, 0l, 1155643113l);
  const yarn3::parameter_type
  yarn3::LEcuyer3=parameter_type(65338l, 0l, 64636l);
  
  // Random number engine concept
  yarn3::yarn3(yarn3::parameter_type P) :
    P(P), S() { }

  yarn3::yarn3(unsigned long s, yarn3::parameter_type P) :
    P(P), S() {
    seed(s);
  }

  void yarn3::seed() {
    (*this)=yarn3();
  }

  void yarn3::seed(unsigned long s) {
    long long t=s;
    t%=modulus;
    if (t<0)
      t+=modulus;
    S.r1=static_cast<result_type>(t);
    S.r2=1;
    S.r3=1;
  }

  void yarn3::seed(yarn3::result_type s1, yarn3::result_type s2,
		   yarn3::result_type s3) {
    S.r1=s1%modulus;
    if (S.r1<0)
      S.r1+=modulus;
    S.r2=s2%modulus;
    if (S.r2<0)
      S.r2+=modulus;
    S.r3=s3%modulus;
    if (S.r3<0)
      S.r3+=modulus;
  }

  // Equality comparable concept
  bool operator==(const yarn3 &R1, const yarn3 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const yarn3 &R1, const yarn3 &R2) {
    return not (R1==R2);
  }

  // Parallel random number generator concept
  void yarn3::split(unsigned int s, unsigned int n) {
    if (s<1 or n>=s)
      throw std::invalid_argument("invalid argument for trng::yarn3::split");
    long q0, q1, q2, q3, q4, q5;
    if (s>1) {
      for (unsigned int i(0); i<=n; ++i) step();
      q0=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q1=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q2=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q3=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q4=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q5=S.r1;
      std::vector<long> a(3), b(9);
      a[0]=q3;  b[0]=q2;  b[1]=q1;  b[2]=q0;
      a[1]=q4;  b[3]=q3;  b[4]=q2;  b[5]=q1;
      a[2]=q5;  b[6]=q4;  b[7]=q3;  b[8]=q2;
      utility::gauss(b, a, modulus);
      P.a1=a[0];  P.a2=a[1];  P.a3=a[2];
      S.r1=q2;    S.r2=q1;    S.r3=q0;
      backward();
      backward();
      backward();
    }
  }
  
  void yarn3::jump2(unsigned int s) {
    std::vector<long> b(9), c(9), d(3), r(3);
    long t1(P.a1), t2(P.a2), t3(P.a3);
    b[0]=P.a1;  b[1]=P.a2;  b[2]=P.a3;
    b[3]=1l;    b[4]=0l;    b[5]=0l;
    b[6]=0l;    b[7]=1l;    b[8]=0l;
    for (unsigned int i(0); i<s; ++i)
      if ((i&1)==0)
	utility::matrix_mult(b, b, c, modulus);
      else
	utility::matrix_mult(c, c, b, modulus);
    r[0]=S.r1;  r[1]=S.r2;  r[2]=S.r3;
    if ((s&1)==0)
      utility::matrix_vec_mult(b, r, d, modulus);
    else
      utility::matrix_vec_mult(c, r, d, modulus);
    S.r1=d[0];  S.r2=d[1];  S.r3=d[2];
    P.a1=t1;    P.a2=t2;    P.a3=t3;
  }

  void yarn3::jump(unsigned long long s) {
    unsigned int i(0);
    while (s>0) {
      if (s%2==1)
	jump2(i);
      ++i;
      s>>=1;
    }
  }

  // other usefull methods
  const char * const yarn3::name_str="yarn3";
  
  const char * yarn3::name() {
    return name_str;
  }

  void yarn3::backward() {
    long t;
    if (P.a3!=0l) {
      t=S.r1;
      t-=static_cast<long>((static_cast<long long>(P.a1)*
			    static_cast<long long>(S.r2))%modulus);
      if (t<0l)
	t+=modulus;
      t-=static_cast<long>((static_cast<long long>(P.a2)*
			    static_cast<long long>(S.r3))%modulus);
    if (t<0l)
      t+=modulus;
    t=static_cast<long>((static_cast<long long>(t)*
                         static_cast<long long>
                         (utility::modulo_invers(P.a3, modulus)))%modulus);
    } else if (P.a2!=0l) {
      t=S.r2;
      t-=static_cast<long>((static_cast<long long>(P.a1)*
			    static_cast<long long>(S.r3))%modulus);
      if (t<0l)
	t+=modulus;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a2, modulus)))%modulus);
    } else if (P.a1!=0l) {
      t=S.r3;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a1, modulus)))%modulus);
    } else
      t=0l;
    S.r1=S.r2;  S.r2=S.r3;  S.r3=t;
  }

  utility::power<yarn3::modulus, yarn3::gen> yarn3::parameter_type::g;

}
