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
#include <trng/yarn2.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/uniform_int_dist.hpp>
#include <trng/exponential_dist.hpp>
#include <trng/normal_dist.hpp>
#include <trng/cauchy_dist.hpp>
#include <trng/logistic_dist.hpp>
#include <trng/lognormal_dist.hpp>
#include <trng/pareto_dist.hpp>
#include <trng/powerlaw_dist.hpp>
#include <trng/tent_dist.hpp>
#include <trng/weibull_dist.hpp>
#include <trng/extreme_value_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <trng/chi_square_dist.hpp>
#include <trng/student_t_dist.hpp>
#include <trng/snedecor_f_dist.hpp>
#include <trng/rayleigh_dist.hpp>
#include <trng/bernoulli_dist.hpp>
#include <trng/binomial_dist.hpp>
#include <trng/geometric_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <trng/discrete_dist.hpp>

template<typename R_type>
void distribution_test() {
  R_type R(1000);
  std::cout << "Generator " << R_type::name() << "\n\n";
  {
    trng::uniform_dist g(50, 100);
    std::cout << "Uniform distribution\n" 
	      << "--------------------\n";
    std::cout << "double   uniform [50, 100)\t" << g(R) << "\n\n";
  }
  {
    trng::uniform01_dist g;
    std::cout << "Uniform01 distribution\n" 
	      << "----------------------\n";
    std::cout << "double   uniform01 [0, 1)\t" << g(R) << "\n\n";
  }
  {
    trng::uniform_int_dist g(98, 100);
    std::cout << "Uniform_int distribution\n" 
	      << "------------------------\n";
    std::cout << "int      uniform_int [98, 100)\t" << g(R) << "\n\n";
  }
  {
    trng::exponential_dist g(5.0);
    std::cout << "Exponential distribution\n" 
	      << "------------------------\n";
    std::cout << "double   exponential (5.0)\t" << g(R) << "\n\n";
  }
  {
    trng::normal_dist g(5.0, 2.0);
    std::cout << "Normal distribution\n" 
	      << "-------------------\n";
    std::cout << "double   normal (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::cauchy_dist g(5.0, 2.0);
    std::cout << "Cauchy distribution\n" 
	      << "-------------------\n";
    std::cout << "double   cauchy (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::logistic_dist g(5.0, 2.0);
    std::cout << "Logistic distribution\n" 
	      << "---------------------\n";
    std::cout << "double   logistic (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::lognormal_dist g(5.0, 2.0);
    std::cout << "Lognormal distribution\n" 
	      << "----------------------\n";
    std::cout << "double   lognormal (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::pareto_dist g(5.0, 2.0);
    std::cout << "Pareto distribution\n" 
	      << "-------------------\n";
    std::cout << "double   pareto (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::powerlaw_dist g(5.0, 2.0);
    std::cout << "Power-law distribution\n" 
	      << "----------------------\n";
    std::cout << "double   powerlaw (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::tent_dist g(5.0, 2.0);
    std::cout << "Tent distribution\n" 
	      << "-----------------\n";
    std::cout << "double   tent (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::weibull_dist g(5.0, 2.0);
    std::cout << "Weibull distribution\n" 
	      << "--------------------\n";
    std::cout << "double   weibull (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::extreme_value_dist g(5.0, 2.0);
    std::cout << "Extreme value distribution\n" 
	      << "--------------------------\n";
    std::cout << "double   extreme_value (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::gamma_dist g(5.0, 2.0);
    std::cout << "Gamma distribution\n" 
	      << "------------------\n";
    std::cout << "double   gamma (5.0, 2.0)\t" << g(R) << "\n\n";
  }
  {
    trng::chi_square_dist g(10);
    std::cout << "Chi-square distribution\n" 
	      << "-----------------------\n";
    std::cout << "double    (10)\t" << g(R) << "\n\n";
  }
  {
    trng::student_t_dist g(10);
    std::cout << "Student-t distribution\n" 
	      << "----------------------\n";
    std::cout << "double    (10)\t" << g(R) << "\n\n";
  }
  {
    trng::snedecor_f_dist g(10, 8);
    std::cout << "F-distribution\n" 
	      << "--------------\n";
    std::cout << "double    (10, 8)\t" << g(R) << "\n\n";
  }
  {
    trng::rayleigh_dist g(4);
    std::cout << "Rayleigh distribution\n" 
	      << "---------------------\n";
    std::cout << "double    (10)\t" << g(R) << "\n\n";
  }
  {
    trng::bernoulli_dist<int> g(0.33333, 1, 2);
    std::cout << "Bernoulli distribution\n" 
	      << "----------------------\n";
    std::cout << "bool    (0.33333)\t" << g(R) << "\n\n";
  }
  {
    trng::binomial_dist g(1.0/3.0, 8);
    std::cout << "Binomial distribution\n" 
	      << "---------------------\n";
    std::cout << "int    (0.33333, 8)\t" << g(R) << "\n\n";
  }
  {
    trng::geometric_dist g(2.0/3.0);
    std::cout << "Geometric distribution\n" 
	      << "----------------------\n";
    std::cout << "int    (0.66667)\t" << g(R) << "\n\n";
  }
  {
    trng::poisson_dist g(10.0);
    std::cout << "Poisson distribution\n" 
	      << "--------------------\n";
    std::cout << "int    (0.66667)\t" << g(R) << "\n\n";
  }
}

int main() {
  distribution_test<trng::yarn2>();
  return EXIT_SUCCESS;
}
