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

#if !(defined TRNG_CONSTANTS_HPP)

#define TRNG_CONSTANTS_HPP

#define TRNG_NEW_CONSTANT(type, value, x) static constexpr type x = value

namespace trng {

  namespace math {

    template<typename T>
    class constants;

    template<>
    class constants<float> {
    public:
      TRNG_NEW_CONSTANT(float, 0.5f, one_half);
      TRNG_NEW_CONSTANT(float, 3.14159265358979323846264f, pi);
      TRNG_NEW_CONSTANT(float, .915965594177219015054604f, catalan);
      TRNG_NEW_CONSTANT(float, 2.71828182845904523536029f, e);
      TRNG_NEW_CONSTANT(float, .577215664901532860606512f, gamma);
      TRNG_NEW_CONSTANT(float, .693147180559945309417232f, ln_2);
      TRNG_NEW_CONSTANT(float, 1.41421356237309504880169f, sqrt_2);
      TRNG_NEW_CONSTANT(float, 2.50662827463100050241577f, sqrt_2pi);
      TRNG_NEW_CONSTANT(float, 0.3183098861837906715377676f, one_over_pi);
      TRNG_NEW_CONSTANT(float, .707106781186547524400845f, one_over_sqrt_2);
      TRNG_NEW_CONSTANT(float, .398942280401432677939946f, one_over_sqrt_2pi);
      TRNG_NEW_CONSTANT(float, .797884560802865355879892f, sqrt_2_over_pi);
      TRNG_NEW_CONSTANT(float, .886226925452758013649085f, sqrt_pi_over_2);
      TRNG_NEW_CONSTANT(float, 1.77245385090551602729817f, sqrt_pi);
    };

    template<>
    class constants<double> {
    public:
      TRNG_NEW_CONSTANT(double, 0.5, one_half);
      TRNG_NEW_CONSTANT(double, 3.14159265358979323846264, pi);
      TRNG_NEW_CONSTANT(double, .915965594177219015054604, catalan);
      TRNG_NEW_CONSTANT(double, 2.71828182845904523536029, e);
      TRNG_NEW_CONSTANT(double, .577215664901532860606512, gamma);
      TRNG_NEW_CONSTANT(double, .693147180559945309417232, ln_2);
      TRNG_NEW_CONSTANT(double, 1.41421356237309504880169, sqrt_2);
      TRNG_NEW_CONSTANT(double, 2.50662827463100050241577, sqrt_2pi);
      TRNG_NEW_CONSTANT(double, 0.3183098861837906715377676, one_over_pi);
      TRNG_NEW_CONSTANT(double, .707106781186547524400845, one_over_sqrt_2);
      TRNG_NEW_CONSTANT(double, .398942280401432677939946, one_over_sqrt_2pi);
      TRNG_NEW_CONSTANT(double, .797884560802865355879892, sqrt_2_over_pi);
      TRNG_NEW_CONSTANT(double, .886226925452758013649085, sqrt_pi_over_2);
      TRNG_NEW_CONSTANT(double, 1.77245385090551602729817, sqrt_pi);
    };

    template<>
    class constants<long double> {
    public:
      TRNG_NEW_CONSTANT(long double, 0.5l, one_half);
      TRNG_NEW_CONSTANT(long double, 3.141592653589793238462643383279502884197l, pi);
      TRNG_NEW_CONSTANT(long double, .9159655941772190150546035149323841107741l, catalan);
      TRNG_NEW_CONSTANT(long double, 2.718281828459045235360287471352662497757l, e);
      TRNG_NEW_CONSTANT(long double, .5772156649015328606065120900824024310422l, gamma);
      TRNG_NEW_CONSTANT(long double, .6931471805599453094172321214581765680755l, ln_2);
      TRNG_NEW_CONSTANT(long double, 1.414213562373095048801688724209698078570l, sqrt_2);
      TRNG_NEW_CONSTANT(long double, 2.506628274631000502415765284811045253008l, sqrt_2pi);
      TRNG_NEW_CONSTANT(long double, .7071067811865475244008443621048490392850l,
                        one_over_sqrt_2);
      TRNG_NEW_CONSTANT(long double, 0.3183098861837906715377675267450287240689l, one_over_pi);
      TRNG_NEW_CONSTANT(long double, .3989422804014326779399460599343818684758l,
                        one_over_sqrt_2pi);
      TRNG_NEW_CONSTANT(long double, .7978845608028653558798921198687637369517l,
                        sqrt_2_over_pi);
      TRNG_NEW_CONSTANT(long double, .8862269254527580136490837416705725913990l,
                        sqrt_pi_over_2);
      TRNG_NEW_CONSTANT(long double, 1.772453850905516027298167483341145182798l, sqrt_pi);
    };

  }  // namespace math

}  // namespace trng

#undef TRNG_NEW_CONSTANT

#endif
