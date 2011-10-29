// Copyright (C) 2000-2007 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#include <trng/mrg5.hpp>

namespace trng {

  // Uniform random number generator concept
  const mrg5::result_type mrg5::modulus=2147483647l;
  const mrg5::result_type mrg5::min=0l;
  const mrg5::result_type mrg5::max=2147483646l;

  // Parameter and status classes

  // EqualityComparable concept
  bool operator==(const mrg5::parameter_type &P1,
                  const mrg5::parameter_type &P2) {
    return P1.a1==P2.a1 && P1.a2==P2.a2 && P1.a3==P2.a3 && P1.a4==P2.a4 && P1.a5==P2.a5;
  }

  bool operator!=(const mrg5::parameter_type &P1,
                  const mrg5::parameter_type &P2) {
    return !(P1==P2);
  }

  // Equality comparable concept
  bool operator==(const mrg5::status_type &S1,
                  const mrg5::status_type &S2) {
    return S1.r1==S2.r1 && S1.r2==S2.r2 && S1.r3==S2.r3 && S1.r4==S2.r4 && S1.r5==S2.r5;
  }

  bool operator!=(const mrg5::status_type &S1,
                  const mrg5::status_type &S2) {
    return !(S1==S2);
  }

  const mrg5::parameter_type
  mrg5::LEcuyer1=parameter_type(107374182l, 0l, 0l, 0l, 104480l);
  
  // Random number engine concept
  mrg5::mrg5(mrg5::parameter_type P) :
    P(P), S() { }

  mrg5::mrg5(unsigned long s, mrg5::parameter_type P) :
    P(P), S() {
    seed(s);
  }

  void mrg5::seed() {
    (*this)=mrg5();
  }

  void mrg5::seed(unsigned long s) {
    long long t=s;
    t%=modulus;
    if (t<0)
      t+=modulus;
    S.r1=static_cast<result_type>(t);
    S.r2=1;
    S.r3=1;
    S.r4=1;
    S.r5=1;
  }

  void mrg5::seed(mrg5::result_type s1, mrg5::result_type s2,
                  mrg5::result_type s3, mrg5::result_type s4, 
		  mrg5::result_type s5) {
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
    S.r5=s5%modulus;
    if (S.r5<0)
      S.r5+=modulus;
  }

  // Equality comparable concept
  bool operator==(const mrg5 &R1, const mrg5 &R2) {
    return R1.P==R2.P && R1.S==R2.S;
  }

  bool operator!=(const mrg5 &R1, const mrg5 &R2) {
    return !(R1==R2);
  }

  // Parallel random number generator concept
  void mrg5::split(unsigned int s, unsigned int n) {
    if (s<1 || n>=s)
      throw std::invalid_argument("invalid argument for trng::mrg5::split");
    long q0, q1, q2, q3, q4, q5, q6, q7, q8, q9;
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
      for (unsigned int i(0); i<s; ++i) step();
      q6=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q7=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q8=S.r1;
      for (unsigned int i(0); i<s; ++i) step();
      q9=S.r1;
      std::vector<long> a(5), b(25);
      a[ 0]=q5;  b[ 0]=q4;  b[ 1]=q3;  b[ 2]=q2;  b[ 3]=q1;  b[ 4]=q0;
      a[ 1]=q6;  b[ 5]=q5;  b[ 6]=q4;  b[ 7]=q3;  b[ 8]=q2;  b[ 9]=q1;
      a[ 2]=q7;  b[10]=q6;  b[11]=q5;  b[12]=q4;  b[13]=q3;  b[14]=q2;
      a[ 3]=q8;  b[15]=q7;  b[16]=q6;  b[17]=q5;  b[18]=q4;  b[19]=q3;
      a[ 4]=q9;  b[20]=q8;  b[21]=q7;  b[22]=q6;  b[23]=q5;  b[24]=q4;
      utility::gauss(b, a, modulus);
      P.a1=a[0];   P.a2=a[1];   P.a3=a[2];   P.a4=a[3];   P.a5=a[4];
      S.r1=q4;     S.r2=q3;     S.r3=q2;     S.r4=q1;     S.r5=q0;
      backward();
      backward();
      backward();
      backward();
      backward();
    }
  }
  
  void mrg5::jump2(unsigned int s) {
    std::vector<long> b(25), c(25), d(5), r(5);
    long t1(P.a1), t2(P.a2), t3(P.a3), t4(P.a4), t5(P.a5);
    b[ 0]=P.a1;  b[ 1]=P.a2;  b[ 2]=P.a3;  b[ 3]=P.a4;  b[ 4]=P.a5;
    b[ 5]=1l;    b[ 6]=0l;    b[ 7]=0l;    b[ 8]=0l;    b[ 9]=0l;
    b[10]=0l;    b[11]=1l;    b[12]=0l;    b[13]=0l;    b[14]=0l;
    b[15]=0l;    b[16]=0l;    b[17]=1l;    b[18]=0l;    b[19]=0l;
    b[20]=0l;    b[21]=0l;    b[22]=0l;    b[23]=1l;    b[24]=0l;
    for (unsigned int i(0); i<s; ++i)
      if ((i&1)==0)
	utility::matrix_mult(b, b, c, modulus);
      else
	utility::matrix_mult(c, c, b, modulus);
    r[0]=S.r1;  r[1]=S.r2;  r[2]=S.r3;  r[3]=S.r4;  r[4]=S.r5;
    if ((s&1)==0)
      utility::matrix_vec_mult(b, r, d, modulus);
    else
      utility::matrix_vec_mult(c, r, d, modulus);
    S.r1=d[0];  S.r2=d[1];  S.r3=d[2];  S.r4=d[3];  S.r5=d[4];
    P.a1=t1;    P.a2=t2;    P.a3=t3;    P.a4=t4;    P.a5=t5;
  }

  void mrg5::jump(unsigned long long s) {
    unsigned int i(0);
    while (s>0) {
      if (s%2==1)
	jump2(i);
      ++i;
      s>>=1;
    }
  }

  // other usefull methods
  const char * const mrg5::name_str="mrg5";
  
  const char * mrg5::name() {
    return name_str;
  }

  void mrg5::backward() {
    long t;
    if (P.a5!=0l) {
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
      t-=static_cast<long>((static_cast<long long>(P.a4)*
                          static_cast<long long>(S.r5))%modulus);
      if (t<0l)
	t+=modulus;
      t=static_cast<long>((static_cast<long long>(t)*
			   static_cast<long long>
			   (utility::modulo_invers(P.a5, modulus)))%modulus);
    } else if (P.a4!=0l) {
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
    S.r1=S.r2;  S.r2=S.r3;  S.r3=S.r4;  S.r4=S.r5;  S.r5=t;
  }
  
}
