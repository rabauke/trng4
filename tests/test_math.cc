// Copyright (c) 2000-2019, Heiko Bauke
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

#include <vector>
#include <ciso646>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/special_functions.hpp>

using floats = boost::mpl::list<float, double, long double>;

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_math)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_special_functions)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_Phi, T, floats) {
  struct test_values {
    const T x, y;
  };
  const std::vector<test_values> values{
      {-8.00000000000000000000000000e+00l, 6.22096057427178412351599517e-16l},
      {-7.00000000000000000000000000e+00l, 1.27981254388583500438362369e-12l},
      {-6.00000000000000000000000000e+00l, 9.86587645037698140700864132e-10l},
      {-5.00000000000000000000000000e+00l, 2.86651571879193911673752333e-07l},
      {-4.00000000000000000000000000e+00l, 3.16712418331199212537707567e-05l},
      {-3.00000000000000000000000000e+00l, 1.34989803163009452665181477e-03l},
      {-2.00000000000000000000000000e+00l, 2.27501319481792072002826372e-02l},
      {-1.00000000000000000000000000e+00l, 1.58655253931457051414767454e-01l},
      {0.00000000000000000000000000e+00l, 5.00000000000000000000000000e-01l},
      {1.00000000000000000000000000e+00l, 8.41344746068542948585232546e-01l},
      {2.00000000000000000000000000e+00l, 9.77249868051820792799717363e-01l},
      {3.00000000000000000000000000e+00l, 9.98650101968369905473348185e-01l},
      {4.00000000000000000000000000e+00l, 9.99968328758166880078746229e-01l},
      {5.00000000000000000000000000e+00l, 9.99999713348428120806088326e-01l},
      {6.00000000000000000000000000e+00l, 9.99999999013412354962301859e-01l},
      {7.00000000000000000000000000e+00l, 9.99999999998720187456114165e-01l},
      {8.00000000000000000000000000e+00l, 9.99999999999999377903942573e-01l}};
  const T eps{32 * std::numeric_limits<T>::epsilon()};
  for (auto &v : values) {
    const T x{v.x};
    T y_min{(1 - eps) * std::abs(v.y)};
    T y_max{(1 + eps) * std::abs(v.y)};
    if (y_min > y_max)
      std::swap(y_min, y_max);
    const T y{trng::math::Phi(x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", Phi(x) = " << y << ", err = " << err
        << ", rel_err = " << rel_err;
    BOOST_TEST((y_min < y and y < y_max), mes.str().c_str());
  }
}
BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
