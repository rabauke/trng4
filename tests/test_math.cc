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
#include <tuple>
#include <utility>
#include <limits>
#include <cmath>
#include <ciso646>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <trng/special_functions.hpp>

#include "type_names.hpp"

using floats = boost::mpl::list<float, double, long double>;

template<typename T>
std::tuple<T, T> bounds(T y) {
  const T eps{32 * std::numeric_limits<T>::epsilon()};
  T y_min{(1 - eps) * std::abs(y)};
  T y_max{(1 + eps) * std::abs(y)};
  if (y_min > y_max)
    std::swap(y_min, y_max);
  if (std::abs(y_min) < 32 * std::numeric_limits<T>::min())
    y_min = -32 * std::numeric_limits<T>::min();
  if (std::abs(y_max) < 32 * std::numeric_limits<T>::min())
    y_max = +32 * std::numeric_limits<T>::min();
  return {y_min, y_max};
}

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_math)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_special_functions)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_asech, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T x, y;
  };
  const std::vector<test_values> values{
      {9.5367431640625000e-07l, 1.45560907917586241240864312e+01l},
      {9.7656250000000000e-04l, 7.62461874774073403685358701e+00l},
      {7.8125000000000000e-03l, 5.54516218534124216700376390e+00l},
      {1.0000000000000000e+00l, 0.00000000000000000000000000e+00l}};
  for (auto &v : values) {
    const T x{v.x};
    const T y{trng::math::asech(x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", asech(x) = " << y << ", err = " << err
        << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_acsch, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T x, y;
  };
  const std::vector<test_values> values{
      {9.5367431640625000e-07l, 1.45560907917590788714373177e+01l},
      {9.7656250000000000e-04l, 7.62461922457789224006893719e+00l},
      {1.0000000000000000e+00l, 8.81373587019543025232609325e-01l},
      {1.0240000000000000e+03l, 9.76562344779637510763910959e-04l},
      {1.0485760000000000e+06l, 9.53674316406105439710335325e-07l}};
  for (auto &v : values) {
    const T x{v.x};
    const T y{trng::math::acsch(x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", acsch(x) = " << y << ", err = " << err
        << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_acoth, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T x, y;
  };
  const std::vector<test_values> values{
      {1.0009765625000000e+00l, 3.81255357411944987738904519e+00l},
      {1.0078125000000000e+00l, 2.77453804244760989917589716e+00l},
      {2.0000000000000000e+00l, 5.49306144334054845697622618e-01l},
      {1.2800000000000000e+02l, 7.81265895154042091032347128e-03l},
      {1.0240000000000000e+03l, 9.76562810441035840964450030e-04l}};
  for (auto &v : values) {
    const T x{v.x};
    const T y{trng::math::acoth(x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", aocth(x) = " << y << ", err = " << err
        << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_special_functions)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_Phi, T, floats) {
  const std::string &name{type_name<T>::name};
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
  for (auto &v : values) {
    const T x{v.x};
    const T y{trng::math::Phi(x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", Phi(x) = " << y << ", err = " << err
        << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_GammaP, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T s, x, y;
  };
  const std::vector<test_values> values{{2.0l, 0.0l, 0.00000000000000000000000000e+00l},
                                        {2.0l, 1.0l, 2.64241117657115356808952460e-01l},
                                        {2.0l, 2.0l, 5.93994150290161924318001515e-01l},
                                        {2.0l, 3.0l, 8.00851726528544228082630337e-01l},
                                        {2.0l, 4.0l, 9.08421805556329098531409894e-01l},
                                        {2.0l, 5.0l, 9.59572318005487197420183709e-01l},
                                        {2.0l, 6.0l, 9.82648734763335491038683828e-01l},
                                        {2.0l, 7.0l, 9.92704944275563870335974911e-01l},
                                        {2.0l, 8.0l, 9.96980836348877393450607498e-01l}};
  for (auto &v : values) {
    const T s{v.s};
    const T x{v.x};
    const T y{trng::math::GammaP(s, x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, s = " << s << ", x = " << x << ", GammaP(s, x) = " << y
        << ", err = " << err << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_GammaQ, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T s, x, y;
  };
  const std::vector<test_values> values{
      {2.0l, 0.0l, 1.00000000000000000000000000e+00l},
      {2.0l, 1.0l, 7.35758882342884643191047540e-01l},
      {2.0l, 2.0l, 4.06005849709838075681998485e-01l},
      {2.0l, 3.0l, 1.99148273471455771917369663e-01l},
      {2.0l, 4.0l, 9.15781944436709014685901064e-02l},
      {2.0l, 5.0l, 4.04276819945128025798162905e-02l},
      {2.0l, 6.0l, 1.73512652366645089613161720e-02l},
      {2.0l, 7.0l, 7.29505572443612966402508868e-03l},
      {2.0l, 8.0l, 3.01916365112260654939250213e-03l},
  };
  for (auto &v : values) {
    const T s{v.s};
    const T x{v.x};
    const T y{trng::math::GammaQ(s, x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, s = " << s << ", x = " << x << ", GammaQ(s, x) = " << y
        << ", err = " << err << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inc_gamma, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T s, x, y;
  };
  const std::vector<test_values> values{{6.0l, 2.0l, 1.98763301767373266843244038e+00l},
                                        {6.0l, 4.0l, 2.57843535563513765691144986e+01l},
                                        {6.0l, 6.0l, 6.65184430362466506643774675e+01l},
                                        {6.0l, 8.0l, 9.70516725504449701299064127e+01l},
                                        {6.0l, 10.0l, 1.11949684454516186125708912e+02l},
                                        {6.0l, 12.0l, 1.17559076469968595452748762e+02l},
                                        {6.0l, 14.0l, 1.19336154036275981201164822e+02l}};
  for (auto &v : values) {
    const T s{v.s};
    const T x{v.x};
    const T y{trng::math::inc_gamma(s, x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, s = " << s << ", x = " << x << ", inc_gamma(s, x) = " << y
        << ", err = " << err << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inc_Gamma, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T s, x, y;
  };
  const std::vector<test_values> values{{6.0l, 2.0l, 1.18012366982326267331567560e+02l},
                                        {6.0l, 4.0l, 9.42156464436486234308855014e+01l},
                                        {6.0l, 6.0l, 5.34815569637533493356225325e+01l},
                                        {6.0l, 8.0l, 2.29483274495550298700935873e+01l},
                                        {6.0l, 10.0l, 8.05031554548381387429108754e+00l},
                                        {6.0l, 12.0l, 2.44092353003140454725123793e+00l},
                                        {6.0l, 14.0l, 6.63845963724018798835178155e-01l}};
  for (auto &v : values) {
    const T s{v.s};
    const T x{v.x};
    const T y{trng::math::inc_Gamma(s, x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, s = " << s << ", x = " << x << ", inc_Gamma(s, x) = " << y
        << ", err = " << err << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_GammaP, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T s, x, y;
  };
  const std::vector<test_values> values{{2.0l, 0.025l, 2.42209278543964902902040257e-01l},
                                        {2.0l, 0.050l, 3.55361510698662052223826033e-01l},
                                        {2.0l, 0.100l, 5.31811608389612020145630298e-01l},
                                        {2.0l, 0.200l, 8.24388309032984609379573008e-01l},
                                        {2.0l, 0.300l, 1.09734921070349164962301800e+00l},
                                        {2.0l, 0.400l, 1.37642134206288679917880939e+00l},
                                        {2.0l, 0.500l, 1.67834699001666065341288451e+00l},
                                        {2.0l, 0.600l, 2.02231324532465691594936386e+00l},
                                        {2.0l, 0.800l, 2.99430834700212208501332324e+00l},
                                        {2.0l, 0.800l, 2.99430834700212208501332324e+00l},
                                        {2.0l, 0.900l, 3.88972016986742905790398022e+00l},
                                        {2.0l, 0.950l, 4.74386451839057837585027379e+00l},
                                        {2.0l, 0.975l, 5.57164339093889859722728276e+00l}};
  for (auto &v : values) {
    const T s{v.s};
    const T x{v.x};
    const T y{trng::math::inv_GammaP(s, x)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, s = " << s << ", x = " << x << ", inv_GammaP(s, x) = " << y
        << ", err = " << err << ", rel_err = " << rel_err << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_Beta_I, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T x, a, b, y;
  };
  const std::vector<test_values> values{
      {0.00000l, 2.0l, 3.0l, 0.00000000000000000000000000e+00l},
      {0.06250l, 2.0l, 3.0l, 2.15301513671875000000000000e-02l},
      {0.12500l, 2.0l, 3.0l, 7.88574218750000000000000000e-02l},
      {0.18750l, 2.0l, 3.0l, 1.61911010742187500000000000e-01l},
      {0.25000l, 2.0l, 3.0l, 2.61718750000000000000000000e-01l},
      {0.31250l, 2.0l, 3.0l, 3.70407104492187500000000000e-01l},
      {0.37500l, 2.0l, 3.0l, 4.81201171875000000000000000e-01l},
      {0.43750l, 2.0l, 3.0l, 5.88424682617187500000000000e-01l},
      {0.50000l, 2.0l, 3.0l, 6.87500000000000000000000000e-01l},
      {0.56250l, 2.0l, 3.0l, 7.74948120117187500000000000e-01l},
      {0.62500l, 2.0l, 3.0l, 8.48388671875000000000000000e-01l},
      {0.68750l, 2.0l, 3.0l, 9.06539916992187500000000000e-01l},
      {0.75000l, 2.0l, 3.0l, 9.49218750000000000000000000e-01l},
      {0.81250l, 2.0l, 3.0l, 9.77340698242187500000000000e-01l},
      {0.87500l, 2.0l, 3.0l, 9.92919921875000000000000000e-01l},
      {0.93750l, 2.0l, 3.0l, 9.99069213867187500000000000e-01l},
      {1.00000l, 2.0l, 3.0l, 1.00000000000000000000000000e+00l}};
  for (auto &v : values) {
    const T x{v.x};
    const T a{v.a};
    const T b{v.b};
    const T y{trng::math::Beta_I(x, a, b)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", a = " << a << ", b = " << b
        << ", Beta_I(x, a, b) = " << y << ", err = " << err << ", rel_err = " << rel_err
        << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_Beta_I, T, floats) {
  const std::string &name{type_name<T>::name};
  struct test_values {
    const T x, a, b, y;
  };
  const std::vector<test_values> values{
      {0.00000l, 2.0l, 3.0l, 0.00000000000000000000000000e+00l},
      {0.06250l, 2.0l, 3.0l, 1.10103992184032790983487819e-01l},
      {0.12500l, 2.0l, 3.0l, 1.61620271098970995979930090e-01l},
      {0.18750l, 2.0l, 3.0l, 2.04338745361132539624626712e-01l},
      {0.25000l, 2.0l, 3.0l, 2.43022083756076302354362773e-01l},
      {0.31250l, 2.0l, 3.0l, 2.79584458771007570188403447e-01l},
      {0.37500l, 2.0l, 3.0l, 3.15090319077922331556821006e-01l},
      {0.43750l, 2.0l, 3.0l, 3.50271211281557842672885757e-01l},
      {0.50000l, 2.0l, 3.0l, 3.85727568132389548275502751e-01l},
      {0.56250l, 2.0l, 3.0l, 4.22038924326916782546717886e-01l},
      {0.62500l, 2.0l, 3.0l, 4.59853597829354432235024343e-01l},
      {0.68750l, 2.0l, 3.0l, 5.00000000000000000000000000e-01l},
      {0.75000l, 2.0l, 3.0l, 5.43678285419080293035191902e-01l},
      {0.81250l, 2.0l, 3.0l, 5.92879386696679638004654292e-01l},
      {0.87500l, 2.0l, 3.0l, 6.51555286733575644057319469e-01l},
      {0.93750l, 2.0l, 3.0l, 7.30452808307633304719173556e-01l},
      {1.00000l, 2.0l, 3.0l, 1.00000000000000000000000000e+00l},
      //
      {0.00000l, 0.25l, 0.2l, 0.00000000000000000000000000e+00l},
      {0.06250l, 0.25l, 0.2l, 3.04837872563118122746607809e-04l},
      {0.12500l, 0.25l, 0.2l, 4.86316680571784400811646727e-03l},
      {0.18750l, 0.25l, 0.2l, 2.43114758223993589384950197e-02l},
      {0.25000l, 0.25l, 0.2l, 7.43127883348396666110669071e-02l},
      {0.31250l, 0.25l, 0.2l, 1.69503878384992273334335178e-01l},
      {0.37500l, 0.25l, 0.2l, 3.13012435978902318709143682e-01l},
      {0.43750l, 0.25l, 0.2l, 4.88092556785491297708294444e-01l},
      {0.50000l, 0.25l, 0.2l, 6.62494795443082887481897957e-01l},
      {0.56250l, 0.25l, 0.2l, 8.05660175039717940471247796e-01l},
      {0.62500l, 0.25l, 0.2l, 9.03590716794438116834527602e-01l}};
  // Maple fails to calculate reference values for larger x
  for (auto &v : values) {
    const T x{v.x};
    const T a{v.a};
    const T b{v.b};
    const T y{trng::math::inv_Beta_I(x, a, b)};
    const T err{std::abs(y - v.y)};
    const T rel_err{err / std::abs(v.y)};
    std::stringstream mes;
    mes << "sufficent accuracy, x = " << x << ", a = " << a << ", b = " << b
        << ", inv_Beta_I(x, a, b) = " << y << ", err = " << err << ", rel_err = " << rel_err
        << " for " << name;
    auto y_min_max = bounds(v.y);
    const T y_min{std::get<0>(y_min_max)};
    const T y_max{std::get<1>(y_min_max)};
    BOOST_TEST((y_min <= y and y <= y_max), mes.str().c_str());
  }
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
