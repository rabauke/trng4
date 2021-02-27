// Copyright (c) 2000-2021, Heiko Bauke
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

#include <array>
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

template<typename T, std::size_t n>
struct args : private std::array<T, n> {
  using base = std::array<T, n>;
  using base::base;
  using base::operator[];
  using base::size;

  args(T x0) : base{x0} { static_assert(n == 1, "wrong number of arguments"); };

  args(T x0, T x1) : base{x0, x1} { static_assert(n == 2, "wrong number of arguments"); };

  args(T x0, T x1, T x2) : base{x0, x1, x2} {
    static_assert(n == 3, "wrong number of arguments");
  };

  template<typename CharT, typename TraitT>
  friend std::basic_ostream<CharT, TraitT> &operator<<(std::basic_ostream<CharT, TraitT> &out,
                                                       const args &x) {
    for (std::size_t i{0}; i < x.size(); ++i) {
      out << x[i];
      if (i + 1 < n)
        out << ", ";
    }
    return out;
  }
};

template<typename T, std::size_t n>
struct arg_res_tuple {
  args<T, n> x;
  T y;
  arg_res_tuple(T x0, T y) : x{x0}, y{y} {
    static_assert(n == 1, "wrong number of arguments");
  };
  arg_res_tuple(T x0, T x1, T y) : x{x0, x1}, y{y} {
    static_assert(n == 2, "wrong number of arguments");
  };
  arg_res_tuple(T x0, T x1, T x2, T y) : x{x0, x1, x2}, y{y} {
    static_assert(n == 3, "wrong number of arguments");
  };
  static constexpr std::size_t size{n};
};

template<typename T>
std::tuple<T, T> bounds(T y) {
  const T tol{std::numeric_limits<T>::digits *
              std::log2(static_cast<T>(std::numeric_limits<T>::radix))};
  const T eps{tol * std::numeric_limits<T>::epsilon()};
  const T min{tol * std::numeric_limits<T>::min()};
  T y_min{(1 - eps) * y};
  T y_max{(1 + eps) * y};
  if (y_min > y_max)
    std::swap(y_min, y_max);
  if (std::abs(y_min) < min)
    y_min = -min;
  if (std::abs(y_max) < min)
    y_max = +min;
  return {y_min, y_max};
}

template<typename T_ref, typename T>
boost::test_tools::predicate_result check_function(const T_ref &x_yref, const T y,
                                                   const char *const function_name) {
  const std::string &name{type_name<T>::name};
  const T y_ref{x_yref.y};
  const T err{std::abs(y - y_ref)};
  const T rel_err{err / std::abs(y_ref)};
  const std::tuple<T, T> y_min_max{bounds(y_ref)};
  const T y_min{std::get<0>(y_min_max)};
  const T y_max{std::get<1>(y_min_max)};
  if (not(y_min <= y and y <= y_max)) {
    boost::test_tools::predicate_result res(false);
    res.message() << "insufficient accuracy, " << function_name << "(" << x_yref.x
                  << ") yields " << y << " with err = " << err << ", rel_err = " << rel_err
                  << " for " << name;
    return res;
  }
  return true;
}

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_math)

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_special_functions)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_asech, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(9.5367431640625000e-07l), T(1.4556090791758624124086431241014046742934e+01l)},
      {T(9.7656250000000000e-04l), T(7.6246187477407340368535870052246141059924e+00l)},
      {T(7.8125000000000000e-03l), T(5.5451621853412421670037638969171713159640e+00l)},
      {T(1.0000000000000000e+00l), T(0.0000000000000000000000000000000000000000e+00l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::asech(x)};
    auto res = check_function(x_yref, y, "asech");
    BOOST_TEST(res);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_acsch, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(9.5367431640625000e-07l), T(1.4556090791759078871437317705133004262544e+01l)},
      {T(9.7656250000000000e-04l), T(7.6246192245778922400689371862651116355740e+00l)},
      {T(1.0000000000000000e+00l), T(8.8137358701954302523260932497979230902816e-01l)},
      {T(1.0240000000000000e+03l), T(9.7656234477963751076391095890851381048162e-04l)},
      {T(1.0485760000000000e+06l), T(9.5367431640610543971033532524003356450374e-07l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::acsch(x)};
    BOOST_TEST(check_function(x_yref, y, "acsch"));
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_acoth, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(1.0009765625000000e+00l), T(3.8125535741194498773890451928662938293404e+00l)},
      {T(1.0078125000000000e+00l), T(2.7745380424476098991758971573806277258074e+00l)},
      {T(2.0000000000000000e+00l), T(5.4930614433405484569762261846126285232375e-01l)},
      {T(1.2800000000000000e+02l), T(7.8126589515404209103234712760401726663588e-03l)},
      {T(1.0240000000000000e+03l), T(9.7656281044103584096445002988532625423842e-04l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::acoth(x)};
    BOOST_TEST(check_function(x_yref, y, "acoth"));
  }
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(test_suite_special_functions)

BOOST_AUTO_TEST_CASE_TEMPLATE(test_Beta, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(5.0e-01l), T(3.0e+00l), T(1.0666666666666666666666666666666666666667e+00l)},
      {T(5.0e-01l), T(6.0e+00l), T(7.3881673881673881673881673881673881673882e-01l)},
      {T(5.0e-01l), T(9.0e+00l), T(5.9907674025321084144613556378262260615202e-01l)},
      {T(5.0e-01l), T(1.2e+01l), T(5.1701948161767787558607319005757557379579e-01l)},
      {T(5.0e-01l), T(1.5e+01l), T(4.6147455340097409155759237684219465774509e-01l)},
      {T(5.0e-01l), T(1.8e+01l), T(4.2068229082552141735623000084148117941069e-01l)},
      {T(5.0e-01l), T(2.1e+01l), T(3.8909005550720098638900842834281307806150e-01l)},
      {T(5.0e-01l), T(2.4e+01l), T(3.6368978435929563418502202954055587596737e-01l)},
      {T(5.0e-01l), T(2.7e+01l), T(3.4269168111048264699306703275016703527243e-01l)},
      {T(5.0e-01l), T(3.0e+01l), T(3.2495542039483022421219156601582621928744e-01l)},
      {T(5.0e-01l), T(3.3e+01l), T(3.0971551267079477083383801932403733576715e-01l)}};
  for (auto &x_yref : values) {
    const T a{x_yref.x[0]};
    const T b{x_yref.x[1]};
    const T y{trng::math::Beta(a, b)};
    BOOST_TEST(check_function(x_yref, y, "Beta"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_GammaP, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(2.0l), T(0.0l), T(0.0000000000000000000000000000000000000000e+00l)},
      {T(2.0l), T(1.0l), T(2.6424111765711535680895245967707826510838e-01l)},
      {T(2.0l), T(2.0l), T(5.9399415029016192431800151508254678977711e-01l)},
      {T(2.0l), T(3.0l), T(8.0085172652854422808263033739975289347320e-01l)},
      {T(2.0l), T(4.0l), T(9.0842180555632909853140989363379378894044e-01l)},
      {T(2.0l), T(5.0l), T(9.5957231800548719742018370946110945450690e-01l)},
      {T(2.0l), T(6.0l), T(9.8264873476333549103868382798428332475945e-01l)},
      {T(2.0l), T(7.0l), T(9.9270494427556387033597491132472573898821e-01l)},
      {T(2.0l), T(8.0l), T(9.9698083634887739345060749786797225082620e-01l)}};
  for (auto &x_yref : values) {
    const T s{x_yref.x[0]};
    const T x{x_yref.x[1]};
    const T y{trng::math::GammaP(s, x)};
    BOOST_TEST(check_function(x_yref, y, "GammaP"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_GammaQ, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(2.0l), T(0.0l), T(1.0000000000000000000000000000000000000000e+00l)},
      {T(2.0l), T(1.0l), T(7.3575888234288464319104754032292173489162e-01l)},
      {T(2.0l), T(2.0l), T(4.0600584970983807568199848491745321022289e-01l)},
      {T(2.0l), T(3.0l), T(1.9914827347145577191736966260024710652680e-01l)},
      {T(2.0l), T(4.0l), T(9.1578194443670901468590106366206211059560e-02l)},
      {T(2.0l), T(5.0l), T(4.0427681994512802579816290538890545493098e-02l)},
      {T(2.0l), T(6.0l), T(1.7351265236664508961316172015716675240545e-02l)},
      {T(2.0l), T(7.0l), T(7.2950557244361296640250886752742610117898e-03l)},
      {T(2.0l), T(8.0l), T(3.0191636511226065493925021320277491737981e-03l)}};
  for (auto &x_yref : values) {
    const T s{x_yref.x[0]};
    const T x{x_yref.x[1]};
    const T y{trng::math::GammaQ(s, x)};
    BOOST_TEST(check_function(x_yref, y, "GammaQ"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inc_gamma, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(6.0l), T(2.0l), T(1.9876330176737326684324403839936002285453e+00l)},
      {T(6.0l), T(4.0l), T(2.5784353556351376569114498570447050061924e+01l)},
      {T(6.0l), T(6.0l), T(6.6518443036246650664377467512699573572856e+01l)},
      {T(6.0l), T(8.0l), T(9.7051672550444970129906412683582859390980e+01l)},
      {T(6.0l), T(10.0l), T(1.1194968445451618612570891246080316579261e+02l)},
      {T(6.0l), T(12.0l), T(1.1755907646996859545274876206518956856766e+02l)},
      {T(6.0l), T(14.0l), T(1.1933615403627598120116482184534387823867e+02l)}};
  for (auto &x_yref : values) {
    const T s{x_yref.x[0]};
    const T x{x_yref.x[1]};
    const T y{trng::math::inc_gamma(s, x)};
    BOOST_TEST(check_function(x_yref, y, "inc_gamma"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inc_Gamma, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(6.0l), T(2.0l), T(1.1801236698232626733156755961600639977145e+02l)},
      {T(6.0l), T(4.0l), T(9.4215646443648623430885501429552949938076e+01l)},
      {T(6.0l), T(6.0l), T(5.3481556963753349335622532487300426427144e+01l)},
      {T(6.0l), T(8.0l), T(2.2948327449555029870093587316417140609020e+01l)},
      {T(6.0l), T(10.0l), T(8.0503155454838138742910875391968342073876e+00l)},
      {T(6.0l), T(12.0l), T(2.4409235300314045472512379348104314323446e+00l)},
      {T(6.0l), T(14.0l), T(6.6384596372401879883517815465612176132626e-01l)},
  };
  for (auto &x_yref : values) {
    const T s{x_yref.x[0]};
    const T x{x_yref.x[1]};
    const T y{trng::math::inc_Gamma(s, x)};
    BOOST_TEST(check_function(x_yref, y, "inc_Gamma"));
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_GammaP, T, floats) {
  const std::vector<arg_res_tuple<T, 2>> values{
      {T(2.0l), T(0.025l), T(2.4220927854396490290204025650292620932214e-01l)},
      {T(2.0l), T(0.050l), T(3.5536151069866205222382603278533247307144e-01l)},
      {T(2.0l), T(0.100l), T(5.3181160838961202014563029774991313268412e-01l)},
      {T(2.0l), T(0.200l), T(8.2438830903298460937957300776559646539963e-01l)},
      {T(2.0l), T(0.300l), T(1.0973492107034916496230180027123214522258e+00l)},
      {T(2.0l), T(0.400l), T(1.3764213420628867991788093855591344225250e+00l)},
      {T(2.0l), T(0.500l), T(1.6783469900166606534128845120945230848245e+00l)},
      {T(2.0l), T(0.600l), T(2.0223132453246569159493638627662623606088e+00l)},
      {T(2.0l), T(0.800l), T(2.9943083470021220850133232356709534186148e+00l)},
      {T(2.0l), T(0.800l), T(2.9943083470021220850133232356709534186148e+00l)},
      {T(2.0l), T(0.900l), T(3.8897201698674290579039802249268070229527e+00l)},
      {T(2.0l), T(0.950l), T(4.7438645183905783758502737858333195550345e+00l)},
      {T(2.0l), T(0.975l), T(5.5716433909388985972272827567407241238352e+00l)}};
  for (auto &x_yref : values) {
    const T s{x_yref.x[0]};
    const T x{x_yref.x[1]};
    const T y{trng::math::inv_GammaP(s, x)};
    BOOST_TEST(check_function(x_yref, y, "inv_GammaP"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_Beta_I, T, floats) {
  const std::vector<arg_res_tuple<T, 3>> values{
      {T(0.00000l), T(2.0l), T(3.0l), T(0.0000000000000000000000000000000000000000e+00l)},
      {T(0.06250l), T(2.0l), T(3.0l), T(2.1530151367187500000000000000000000000000e-02l)},
      {T(0.12500l), T(2.0l), T(3.0l), T(7.8857421875000000000000000000000000000000e-02l)},
      {T(0.18750l), T(2.0l), T(3.0l), T(1.6191101074218750000000000000000000000000e-01l)},
      {T(0.25000l), T(2.0l), T(3.0l), T(2.6171875000000000000000000000000000000000e-01l)},
      {T(0.31250l), T(2.0l), T(3.0l), T(3.7040710449218750000000000000000000000000e-01l)},
      {T(0.37500l), T(2.0l), T(3.0l), T(4.8120117187500000000000000000000000000000e-01l)},
      {T(0.43750l), T(2.0l), T(3.0l), T(5.8842468261718750000000000000000000000000e-01l)},
      {T(0.50000l), T(2.0l), T(3.0l), T(6.8750000000000000000000000000000000000000e-01l)},
      {T(0.56250l), T(2.0l), T(3.0l), T(7.7494812011718750000000000000000000000000e-01l)},
      {T(0.62500l), T(2.0l), T(3.0l), T(8.4838867187500000000000000000000000000000e-01l)},
      {T(0.68750l), T(2.0l), T(3.0l), T(9.0653991699218750000000000000000000000000e-01l)},
      {T(0.75000l), T(2.0l), T(3.0l), T(9.4921875000000000000000000000000000000000e-01l)},
      {T(0.81250l), T(2.0l), T(3.0l), T(9.7734069824218750000000000000000000000000e-01l)},
      {T(0.87500l), T(2.0l), T(3.0l), T(9.9291992187500000000000000000000000000000e-01l)},
      {T(0.93750l), T(2.0l), T(3.0l), T(9.9906921386718750000000000000000000000000e-01l)},
      {T(1.00000l), T(2.0l), T(3.0l), T(1.0000000000000000000000000000000000000000e+00l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T a{x_yref.x[1]};
    const T b{x_yref.x[2]};
    const T y{trng::math::Beta_I(x, a, b)};
    BOOST_TEST(check_function(x_yref, y, "Beta_I"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_Beta_I, T, floats) {
  const std::vector<arg_res_tuple<T, 3>> values{
      {T(0.00000l), T(2.00l), T(3.0l), T(0.0000000000000000000000000000000000000000e+00l)},
      {T(0.06250l), T(2.00l), T(3.0l), T(1.1010399218403279098348781939518564840313e-01l)},
      {T(0.12500l), T(2.00l), T(3.0l), T(1.6162027109897099597993009042954237965715e-01l)},
      {T(0.18750l), T(2.00l), T(3.0l), T(2.0433874536113253962462671167384161533153e-01l)},
      {T(0.25000l), T(2.00l), T(3.0l), T(2.4302208375607630235436277345892655396111e-01l)},
      {T(0.31250l), T(2.00l), T(3.0l), T(2.7958445877100757018840344667317414747367e-01l)},
      {T(0.37500l), T(2.00l), T(3.0l), T(3.1509031907792233155682100586881160006454e-01l)},
      {T(0.43750l), T(2.00l), T(3.0l), T(3.5027121128155784267288575728522599259027e-01l)},
      {T(0.50000l), T(2.00l), T(3.0l), T(3.8572756813238954827550275115735326717790e-01l)},
      {T(0.56250l), T(2.00l), T(3.0l), T(4.2203892432691678254671788568232128206958e-01l)},
      {T(0.62500l), T(2.00l), T(3.0l), T(4.5985359782935443223502434338992941820028e-01l)},
      {T(0.68750l), T(2.00l), T(3.0l), T(5.0000000000000000000000000000000000000000e-01l)},
      {T(0.75000l), T(2.00l), T(3.0l), T(5.4367828541908029303519190199561893147434e-01l)},
      {T(0.81250l), T(2.00l), T(3.0l), T(5.9287938669667963800465429248808384315208e-01l)},
      {T(0.87500l), T(2.00l), T(3.0l), T(6.5155528673357564405731946917556722675650e-01l)},
      {T(0.93750l), T(2.00l), T(3.0l), T(7.3045280830763330471917355600024013149773e-01l)},
      {T(1.00000l), T(2.00l), T(3.0l), T(1.0000000000000000000000000000000000000000e+00l)},
      //
      {T(0.00000l), T(0.25l), T(0.2l), T(0.0000000000000000000000000000000000000000e+00l)},
      {T(0.06250l), T(0.25l), T(0.2l), T(3.0483787256311812274660780866204614847383e-04l)},
      {T(0.12500l), T(0.25l), T(0.2l), T(4.8631668057178440081164672713718345527002e-03l)},
      {T(0.18750l), T(0.25l), T(0.2l), T(2.4311475822399358938495019664518017293964e-02l)},
      {T(0.25000l), T(0.25l), T(0.2l), T(7.4312788334839666611066907032802850847530e-02l)},
      {T(0.31250l), T(0.25l), T(0.2l), T(1.6950387838499227333433517808491955383025e-01l)},
      {T(0.37500l), T(0.25l), T(0.2l), T(3.1301243597890231870914368223653948288960e-01l)},
      {T(0.43750l), T(0.25l), T(0.2l), T(4.8809255678549129770829444431575692406349e-01l)},
      {T(0.50000l), T(0.25l), T(0.2l), T(6.6249479544308288748189795699313604267852e-01l)},
      {T(0.56250l), T(0.25l), T(0.2l), T(8.0566017503971794047124779538707465925978e-01l)},
      {T(0.62500l), T(0.25l), T(0.2l), T(9.0359071679443811683452760211535163879284e-01l)}};
  // Maple fails to calculate reference values for larger x
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T a{x_yref.x[1]};
    const T b{x_yref.x[2]};
    const T y{trng::math::inv_Beta_I(x, a, b)};
    BOOST_TEST(check_function(x_yref, y, "inv_Beta_I"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_Phi, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(-8.0e+00l), T(6.2209605742717841235159951725881884224887e-16l)},
      {T(-7.0e+00l), T(1.2798125438858350043836236907808329980328e-12l)},
      {T(-6.0e+00l), T(9.8658764503769814070086413239804201866979e-10l)},
      {T(-5.0e+00l), T(2.8665157187919391167375233287464535385442e-07l)},
      {T(-4.0e+00l), T(3.1671241833119921253770756722151298443833e-05l)},
      {T(-3.0e+00l), T(1.3498980316300945266518147675949773778294e-03l)},
      {T(-2.0e+00l), T(2.2750131948179207200282637166533437471776e-02l)},
      {T(-1.0e+00l), T(1.5865525393145705141476745436796207752209e-01l)},
      {T(0.0e+00l), T(5.0000000000000000000000000000000000000000e-01l)},
      {T(1.0e+00l), T(8.4134474606854294858523254563203792247791e-01l)},
      {T(2.0e+00l), T(9.7724986805182079279971736283346656252822e-01l)},
      {T(3.0e+00l), T(9.9865010196836990547334818523240502262217e-01l)},
      {T(4.0e+00l), T(9.9996832875816688007874622924327784870156e-01l)},
      {T(5.0e+00l), T(9.9999971334842812080608832624766712535465e-01l)},
      {T(6.0e+00l), T(9.9999999901341235496230185929913586760196e-01l)},
      {T(7.0e+00l), T(9.9999999999872018745611416499561637630922e-01l)},
      {T(8.0e+00l), T(9.9999999999999937790394257282158764840048e-01l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::Phi(x)};
    BOOST_TEST(check_function(x_yref, y, "Phi"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_Phi, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(9.765625000e-04l), T(-3.0972690781987844623648304970552534107624e+00l)},
      {T(6.2500e-02l), T(-1.5341205443525463117083990590371655257643e+00l)},
      {T(1.2500e-01l), T(-1.1503493803760081782967653108305853365444e+00l)},
      {T(1.8750e-01l), T(-8.8714655901887605668566425981021647661160e-01l)},
      {T(2.5000e-01l), T(-6.7448975019608174320222701454130718538690e-01l)},
      {T(3.1250e-01l), T(-4.8877641111466949891088578294905527635290e-01l)},
      {T(3.7500e-01l), T(-3.1863936396437516302194846367007464334964e-01l)},
      {T(4.3750e-01l), T(-1.5731068461017069552237071807629031362039e-01l)},
      {T(5.0000e-01l), T(0.0000000000000000000000000000000000000000e+00l)},
      {T(5.6250e-01l), T(1.5731068461017069552237071807629031362039e-01l)},
      {T(6.2500e-01l), T(3.1863936396437516302194846367007464334964e-01l)},
      {T(6.8750e-01l), T(4.8877641111466949891088578294905527635290e-01l)},
      {T(7.5000e-01l), T(6.7448975019608174320222701454130718538690e-01l)},
      {T(8.1250e-01l), T(8.8714655901887605668566425981021647661160e-01l)},
      {T(8.7500e-01l), T(1.1503493803760081782967653108305853365444e+00l)},
      {T(9.3750e-01l), T(1.5341205443525463117083990590371655257643e+00l)},
      {T(9.990234375e-01l), T(3.0972690781987844623648304970552534107624e+00l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::inv_Phi(x)};
    BOOST_TEST(check_function(x_yref, y, "inv_Phi(x)"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_erf, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(-9.990234375e-01l), T(-2.3314677736219476723205459143501501091122e+00l)},
      {T(-9.3750e-01l), T(-1.3171503349861307488839297920844487996026e+00l)},
      {T(-8.1250e-01l), T(-9.3197444316109708936320065864204514744200e-01l)},
      {T(-6.8750e-01l), T(-7.1417089760812834681405311779778977078866e-01l)},
      {T(-5.6250e-01l), T(-5.4901309236850152043290414369110772114039e-01l)},
      {T(-4.3750e-01l), T(-4.0950827913413152004135823206066332461754e-01l)},
      {T(-3.1250e-01l), T(-2.8443374892172365391933605508574178560997e-01l)},
      {T(-1.8750e-01l), T(-1.6772722001813860117802187451910453017159e-01l)},
      {T(-6.2500e-02l), T(-5.5445948772782020298989375535954031087215e-02l)},
      {T(6.2500e-02l), T(5.5445948772782020298989375535954031087215e-02l)},
      {T(1.8750e-01l), T(1.6772722001813860117802187451910453017159e-01l)},
      {T(3.1250e-01l), T(2.8443374892172365391933605508574178560997e-01l)},
      {T(4.3750e-01l), T(4.0950827913413152004135823206066332461754e-01l)},
      {T(5.6250e-01l), T(5.4901309236850152043290414369110772114039e-01l)},
      {T(6.8750e-01l), T(7.1417089760812834681405311779778977078866e-01l)},
      {T(8.1250e-01l), T(9.3197444316109708936320065864204514744200e-01l)},
      {T(9.3750e-01l), T(1.3171503349861307488839297920844487996026e+00l)},
      {T(9.990234375e-01l), T(2.3314677736219476723205459143501501091122e+00l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::inv_erf(x)};
    BOOST_TEST(check_function(x_yref, y, "inv_erf"));
  }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(test_inv_erfc, T, floats) {
  const std::vector<arg_res_tuple<T, 1>> values{
      {T(9.765625000000e-04l), T(2.3314677736219476723205459143501501091122e+00l)},
      {T(6.2500e-02l), T(1.3171503349861307488839297920844487996026e+00l)},
      {T(1.8750e-01l), T(9.3197444316109708936320065864204514744200e-01l)},
      {T(3.1250e-01l), T(7.1417089760812834681405311779778977078866e-01l)},
      {T(4.3750e-01l), T(5.4901309236850152043290414369110772114039e-01l)},
      {T(5.6250e-01l), T(4.0950827913413152004135823206066332461754e-01l)},
      {T(6.8750e-01l), T(2.8443374892172365391933605508574178560997e-01l)},
      {T(8.1250e-01l), T(1.6772722001813860117802187451910453017159e-01l)},
      {T(9.3750e-01l), T(5.5445948772782020298989375535954031087215e-02l)},
      {T(1.0625e+00l), T(-5.5445948772782020298989375535954031087215e-02l)},
      {T(1.1875e+00l), T(-1.6772722001813860117802187451910453017159e-01l)},
      {T(1.3125e+00l), T(-2.8443374892172365391933605508574178560997e-01l)},
      {T(1.4375e+00l), T(-4.0950827913413152004135823206066332461754e-01l)},
      {T(1.5625e+00l), T(-5.4901309236850152043290414369110772114039e-01l)},
      {T(1.6875e+00l), T(-7.1417089760812834681405311779778977078866e-01l)},
      {T(1.8125e+00l), T(-9.3197444316109708936320065864204514744200e-01l)},
      {T(1.9375e+00l), T(-1.3171503349861307488839297920844487996026e+00l)},
      {T(1.999023437500e+00l), T(-2.3314677736219476723205459143501501091122e+00l)}};
  for (auto &x_yref : values) {
    const T x{x_yref.x[0]};
    const T y{trng::math::inv_erfc(x)};
    BOOST_TEST(check_function(x_yref, y, "inv_erfc"));
  }
}

BOOST_AUTO_TEST_SUITE_END()

//-----------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
