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

#include <trng/utility.hpp>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>

namespace trng {

  namespace utility {
    
    long modulo_invers(long a, long m) {
      if (a<=0l or m<=1)
	throw std::invalid_argument("invalid argument in trng::utility::modulo_invers");
      long temp, q, flast(0), f(1), m1(m);
      while (a>1) {
	temp=m1%a;
	q=m1/a;
	m1=a;  a=temp;  temp=f;
	f=flast-q*f;
	flast=temp;
      }
      if (a==0)
	throw std::runtime_error("no inversive in trng::utility::modulo_invers");
      return f>=0 ? f : f+m;
    }
    
    //------------------------------------------------------------------

    void gauss(std::vector<long> &a, std::vector<long> &b, long m) {
      if (a.size()!=b.size()*b.size() or a.size()==0 or b.size()==0)
	throw std::invalid_argument("wrong matrix size in trng::utility::gauss");
      // initialize indices
      int n(b.size());
      int rank(0);
      std::vector<long> p(n);
      for (int i(0); i<n; ++i)
	p[i]=i;
      // make matrix triangular
      for (int i(0); i<n; ++i) {
	// search for a pivot element
	if (a[n*p[i]+i]==0l) {
	  // swap rows
	  int j=i+1;
	  while (j<n and a[n*p[j]+i]==0l)
	    j++;
	  if (j<n) {
	    long t=p[i];  p[i]=p[j];  p[j]=t;
	  }
	}
	// is rank small?
	if (a[n*p[i]+i]==0l)
	  break;
	++rank;
	long t=modulo_invers(a[n*p[i]+i], m);
	for (int j=i; j<n; ++j)
	  a[n*p[i]+j]=static_cast<long>
	    ((static_cast<long long>(a[n*p[i]+j])*
	      static_cast<long long>(t))%m);
	b[p[i]]=static_cast<long>
	  ((static_cast<long long>(b[p[i]])*
	    static_cast<long long>(t))%m);
	for (int j=i+1; j<n; ++j) {
	  if (a[n*p[j]+i]!=0l) {
	    t=modulo_invers(a[n*p[j]+i], m);
	    for (int k=i; k<n; ++k) {
	      a[n*p[j]+k]=
		static_cast<long>
		((static_cast<long long>(a[n*p[j]+k])*
		  static_cast<long long>(t))%m);
	      a[n*p[j]+k]-=a[n*p[i]+k];
	      if (a[n*p[j]+k]<0l)
		a[n*p[j]+k]+=m;
	    }
	    b[p[j]]=static_cast<long>
	      ((static_cast<long long>(b[p[j]])*
		static_cast<long long>(t))%m);
	    b[p[j]]-=b[p[i]];
	    if (b[p[j]]<0l)
	      b[p[j]]+=m;
	  }
	}
      }
      // test if a solution exists
      for (int i=rank; i<n; ++i)
	if (b[p[i]]!=0l)
	  throw std::runtime_error("equations system has no solution trng::utility::gauss");
      // solve triangular system
      for (int i=n-2; i>=0; --i)
	for (int j=i+1; j<n; ++j) {
	  b[p[i]]-=static_cast<long>
	    ((static_cast<long long>(a[n*p[i]+j])*
	      static_cast<long long>(b[p[j]]))%m);
	  if (b[p[i]]<0l)
	    b[p[i]]+=m;
	}
      // sort
      for (int i(0); i<n; ++i)
	p[i]=b[p[i]];
      for (int i(0); i<n; ++i)
	b[i]=p[i];
    }
    
    //------------------------------------------------------------------

    void matrix_mult(const std::vector<long> &a, const std::vector<long> &b,
		     std::vector<long> &c, long m) {
      if (a.size()!=b.size())
	throw std::invalid_argument("different sized matrices in trng::utility::matrix_mult");
      int n(static_cast<int>(std::sqrt(static_cast<double>(b.size()))));
      c.resize(n*n);
      for (int i(0); i<n; ++i)
	for (int j(0); j<n; ++j) {
	  long long t(0ll);
	  for (int k(0); k<n; ++k) {
	    t+=(static_cast<long long>(a[j*n+k])*
		static_cast<long long>(b[k*n+i]))%m;
	    if (t>=m)
	      t-=m;
	  }
	  c[j*n+i]=static_cast<long>(t);
	}
    }

    //------------------------------------------------------------------

    void matrix_vec_mult(const std::vector<long> &a, const std::vector<long> &b,
			 std::vector<long> &c, long m) {
      if (a.size()!=b.size()*b.size())
	throw std::invalid_argument("different sized vectors in trng::utility::matrix_vec_mult");
      int n(b.size());
      c.resize(n);
      for (int j(0); j<n; ++j) {
	long long t(0ll);
	for (int k(0); k<n; ++k) {
	  t+=(static_cast<long long>(a[j*n+k])*
	      static_cast<long long>(b[k]))%m;
	  if (t>=m)
	    t-=m;
	}
	c[j]=static_cast<long>(t);
      }
    }
        
  }

}
