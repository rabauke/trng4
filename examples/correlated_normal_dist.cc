// Copyright (C) 2010 Heiko Bauke <heiko.bauke@physics.ox.ac.uk>
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
#include <iomanip>
#include <vector>
#include <trng/lcg64.hpp>
#include <trng/correlated_normal_dist.hpp>

double covariance(const std::vector<double> &v1, const std::vector<double> &v2) {
  std::vector<double>::size_type n=v1.size();
  double m1=0.0, m2=0.0, c=0.0;
  for (std::vector<double>::size_type i=0; i<n; ++i) {
    m1+=v1[i]/n;  m2+=v2[i]/n;
  }
  for (std::vector<double>::size_type i=0; i<n; ++i) 
    c+=(v1[i]-m1)*(v2[i]-m2)/n;
  return c;
}

int main() {
  const int d=4;
  // covariance matrix
  double sig[d][d] = { { 2.0, -0.5,  0.3, -0.3},
		       {-0.5,  3.0, -0.3,  0.3},
		       { 0.3, -0.3,  1.0, -0.3},
		       {-0.3,  0.3, -0.3,  1.0} };
  trng::correlated_normal_dist<> D(&sig[0][0], &sig[d-1][d-1]+1);
  trng::lcg64 R;

  std::vector<double> x1, x2, x3, x4;
  // generate 4-tuples of correlated normal variables
  for (int i=0; i<1000000; ++i) {
    x1.push_back(D(R));  
    x2.push_back(D(R));  
    x3.push_back(D(R)); 
    x4.push_back(D(R));
  }
  // print (empirical) covariance matrix
  std::cout << std::setprecision(4)
	    << covariance(x1, x1) << '\t' << covariance(x1, x2) << '\t'
	    << covariance(x1, x3) << '\t' << covariance(x1, x4) << '\n'
	    << covariance(x2, x1) << '\t' << covariance(x2, x2) << '\t'
	    << covariance(x2, x3) << '\t' << covariance(x2, x4) << '\n'
	    << covariance(x3, x1) << '\t' << covariance(x3, x2) << '\t'
	    << covariance(x3, x3) << '\t' << covariance(x3, x4) << '\n'
	    << covariance(x4, x1) << '\t' << covariance(x4, x2) << '\t'
	    << covariance(x4, x3) << '\t' << covariance(x4, x4) << '\n';
  return EXIT_SUCCESS;
}
