// Copyright (C) 2000-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#include <trng/mrg4.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // EqualityComparable concept
  bool operator==(const mrg4::parameter_type &P1,
                  const mrg4::parameter_type &P2) {
    return P1.a1==P2.a1 and P1.a2==P2.a2 and P1.a3==P2.a3 and P1.a4==P2.a4;
  }

  bool operator!=(const mrg4::parameter_type &P1,
                  const mrg4::parameter_type &P2) {
    return not (P1==P2);
  }

  // Equality comparable concept
  bool operator==(const mrg4::status_type &S1,
                  const mrg4::status_type &S2) {
    return S1.r1==S2.r1 and S1.r2==S2.r2 and S1.r3==S2.r3 and S1.r4==S2.r4;
  }

  bool operator!=(const mrg4::status_type &S1,
                  const mrg4::status_type &S2) {
    return not (S1==S2);
  }

  const mrg4::parameter_type
  mrg4::LEcuyer1=parameter_type(2001982722l, 1412284257l, 1155380217l, 1668339922l);
  const mrg4::parameter_type
  mrg4::LEcuyer2=parameter_type(64886l, 0l, 0l, 64322l);
  
  // Random number engine concept
  mrg4::mrg4(mrg4::parameter_type P) :
    P(P), S() { }

  mrg4::mrg4(unsigned long s, mrg4::parameter_type P) :
    P(P), S() {
    seed(s);
  }

  void mrg4::seed() {
    (*this)=mrg4();
  }

  void mrg4::seed(unsigned long s) {
    long long t=s;
    t%=modulus;
    if (t<0)
      t+=modulus;
    S.r1=static_cast<result_type>(t);
    S.r2=1;
    S.r3=1;
    S.r4=1;
  }

  void mrg4::seed(mrg4::result_type s1, mrg4::result_type s2,
                  mrg4::result_type s3, mrg4::result_type s4) {
    S.r1=s1%modulus;
    if (S.r1<0)
      S.r1+=modulus;
    S.r2=s2%modulus;
    if (S.r2<0)
      S.r2+=modulus;
    S.r3=s3%modulus;
    if (S.r3<0)
      S.r3+=modulus;
    S.r4=s4%modulus;
    if (S.r4<0)
      S.r4+=modulus;
  }

  // Equality comparable concept
  bool operator==(const mrg4 &R1, const mrg4 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const mrg4 &R1, const mrg4 &R2) {
    return not (R1==R2);
  }

  // Parallel random number generator concept
  void mrg4::split(unsigned int s, unsigned int n) {
    if (s<1 or n>=s)
      throw std::invalid_argument("invalid argument for trng::mrg4::split");
    if (s>1) {
      jump(n+1);  long q0=S.r1;
      jump(s);    long q1=S.r1;
      jump(s);    long q2=S.r1;
      jump(s);    long q3=S.r1;
      jump(s);    long q4=S.r1;
      jump(s);    long q5=S.r1;
      jump(s);    long q6=S.r1;
      jump(s);    long q7=S.r1;
      std::vector<long> a(4), b(16);
      a[ 0]=q4;  b[ 0]=q3;  b[ 1]=q2;  b[ 2]=q1;  b[ 3]=q0;
      a[ 1]=q5;  b[ 4]=q4;  b[ 5]=q3;  b[ 6]=q2;  b[ 7]=q1;
      a[ 2]=q6;  b[ 8]=q5;  b[ 9]=q4;  b[10]=q3;  b[11]=q2;
      a[ 3]=q7;  b[12]=q6;  b[13]=q5;  b[14]=q4;  b[15]=q3;
      utility::gauss(b, a, modulus);
      P.a1=a[0];   P.a2=a[1];   P.a3=a[2];   P.a4=a[3];
      S.r1=q3;     S.r2=q2;     S.r3=q1;     S.r4=q0;
      backward();
      backward();
      backward();
      backward();
    }
  }
  
  void mrg4::jump2(unsigned int s) {
    std::vector<long> b(16), c(16), d(4), r(4);
    long t1(P.a1), t2(P.a2), t3(P.a3), t4(P.a4);
    b[ 0]=P.a1;  b[ 1]=P.a2;  b[ 2]=P.a3;  b[ 3]=P.a4;
    b[ 4]=1l;    b[ 5]=0l;    b[ 6]=0l;    b[ 7]=0l;
    b[ 8]=0l;    b[ 9]=1l;    b[10]=0l;    b[11]=0l;
    b[12]=0l;    b[13]=0l;    b[14]=1l;    b[15]=0l;
    for (unsigned int i(0); i<s; ++i)
      if ((i&1)==0)
	utility::matrix_mult(b, b, c, modulus);
    else
      utility::matrix_mult(c, c, b, modulus);
    r[0]=S.r1;  r[1]=S.r2;  r[2]=S.r3;  r[3]=S.r4;
    if ((s&1)==0)
      utility::matrix_vec_mult(b, r, d, modulus);
    else
      utility::matrix_vec_mult(c, r, d, modulus);
    S.r1=d[0];  S.r2=d[1];  S.r3=d[2];  S.r4=d[3];
    P.a1=t1;    P.a2=t2;    P.a3=t3;    P.a4=t4;
  }

  void mrg4::jump(unsigned long long s) {
    if (s<16) {
      for (unsigned int i(0); i<s; ++i) 
	step();
    } else {
      unsigned int i(0);
      while (s>0) {
	if (s%2==1)
	  jump2(i);
	++i;
	s>>=1;
      }
    }
  }

  // other useful methods
  const char * const mrg4::name_str="mrg4";
  
  const char * mrg4::name() {
    return name_str;
  }

  void mrg4::backward() {
    long t;
    if (P.a4!=0l) {
      t=S.r1;
      t-=static_cast<long>((static_cast<long long>(P.a1)*
			    static_cast<long long>(S.r2))%modulus);
      if (t<0l)
	t+=modulus;
      t-=static_cast<long>((static_cast<long long>(P.a2)*
			    static_cast<long long>(S.r3))%modulus);
      if (t<0l)
	t+=modulus;
      t-=static_cast<long>((static_cast<long long>(P.a3)*
                          static_cast<long long>(S.r4))%modulus);
      if (t<0l)
	t+=modulus;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a4, modulus)))%modulus);
    } else if (P.a3!=0l) {
      t=S.r2;
      t-=static_cast<long>((static_cast<long long>(P.a1)*
			    static_cast<long long>(S.r3))%modulus);
      if (t<0l)
	t+=modulus;
      t-=static_cast<long>((static_cast<long long>(P.a2)*
			    static_cast<long long>(S.r4))%modulus);
      if (t<0l)
	t+=modulus;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a3, modulus)))%modulus);
    } else if (P.a2!=0l) {
      t=S.r3;
      t-=static_cast<long>((static_cast<long long>(P.a1)*
			    static_cast<long long>(S.r4))%modulus);
      if (t<0l)
	t+=modulus;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a2, modulus)))%modulus);
    } else if (P.a1!=0l) {
      t=S.r4;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a1, modulus)))%modulus);
    } else
      t=0l;
    S.r1=S.r2;  S.r2=S.r3;  S.r3=S.r4;  S.r4=t;
  }
  
}
