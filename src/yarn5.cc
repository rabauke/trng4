// Copyright (c) 2000-2010, Heiko Bauke
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.  
// 
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.  
// 
//   * Neither the name of the copyright holder nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.

#include <trng/yarn5.hpp>

namespace trng {

  // Uniform random number generator concept

  // Parameter and status classes

  // EqualityComparable concept
  bool operator==(const yarn5::parameter_type &P1,
                  const yarn5::parameter_type &P2) {
    return P1.a1==P2.a1 and P1.a2==P2.a2 and P1.a3==P2.a3 and P1.a4==P2.a4 and P1.a5==P2.a5;
  }

  bool operator!=(const yarn5::parameter_type &P1,
                  const yarn5::parameter_type &P2) {
    return not (P1==P2);
  }

  // Equality comparable concept
  bool operator==(const yarn5::status_type &S1,
                  const yarn5::status_type &S2) {
    return S1.r1==S2.r1 and S1.r2==S2.r2 and S1.r3==S2.r3 and S1.r4==S2.r4 and S1.r5==S2.r5;
  }

  bool operator!=(const yarn5::status_type &S1,
                  const yarn5::status_type &S2) {
    return not (S1==S2);
  }
  
  const yarn5::parameter_type
  yarn5::LEcuyer1=parameter_type(107374182l, 0l, 0l, 0l, 104480l);

  // Random number engine concept
  yarn5::yarn5(yarn5::parameter_type P) :
    P(P), S() { }

  yarn5::yarn5(unsigned long s, yarn5::parameter_type P) :
    P(P), S() {
    seed(s);
  }

  void yarn5::seed() {
    (*this)=yarn5();
  }

  void yarn5::seed(unsigned long s) {
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

  void yarn5::seed(yarn5::result_type s1, yarn5::result_type s2,
                   yarn5::result_type s3, yarn5::result_type s4, 
		   yarn5::result_type s5) {
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
  bool operator==(const yarn5 &R1, const yarn5 &R2) {
    return R1.P==R2.P and R1.S==R2.S;
  }

  bool operator!=(const yarn5 &R1, const yarn5 &R2) {
    return not (R1==R2);
  }

  // Parallel random number generator concept
  void yarn5::split(unsigned int s, unsigned int n) {
    if (s<1 or n>=s)
      throw std::invalid_argument("invalid argument for trng::yarn5::split");
    if (s>1) {
      jump(n+1);  long q0=S.r1;
      jump(s);    long q1=S.r1;
      jump(s);    long q2=S.r1;
      jump(s);    long q3=S.r1;
      jump(s);    long q4=S.r1;
      jump(s);    long q5=S.r1;
      jump(s);    long q6=S.r1;
      jump(s);    long q7=S.r1;
      jump(s);    long q8=S.r1;
      jump(s);    long q9=S.r1;
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
  
  void yarn5::jump2(unsigned int s) {
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

  void yarn5::jump(unsigned long long s) {
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
  const char * const yarn5::name_str="yarn5";
  
  const char * yarn5::name() {
    return name_str;
  }

  void yarn5::backward() {
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

  utility::power<yarn5::modulus, yarn5::gen> yarn5::parameter_type::g;

}
