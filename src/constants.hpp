#if !(defined TRNG_CONSTANTS_HPP)

#define TRNG_CONSTANTS_HPP

#define TRNG_NEW_CONSTANT(type, value, x) \
    static type x() throw () { \
      return value; \
    }

namespace trng {

  namespace math {  

    template<typename T>
    class constants {
    public:
      static T pi() throw () {
	return T(0);
      }
      static T catalan() throw () {
	return T(0);
      }
      static T e() throw () {
	return T(0);
      }
      static T gamma() throw () {
	return T(0);
      }
      static T ln_2() throw () {
	return T(0);
      }
      static T sqrt_2() throw () {
	return T(0);
      }
      static T one_over_sqrt_2() throw () {
	return T(0);
      }
      static T one_over_sqrt_2pi() throw () {
	return T(0);
      }
    };

    template<> class constants<float> {
    public:
      TRNG_NEW_CONSTANT(float, 3.14159265358979323846264f, pi);
      TRNG_NEW_CONSTANT(float, .915965594177219015054604f, catalan);
      TRNG_NEW_CONSTANT(float, 2.71828182845904523536029f, e);
      TRNG_NEW_CONSTANT(float, .577215664901532860606512f, gamma);
      TRNG_NEW_CONSTANT(float, .693147180559945309417232f, ln_2);
      TRNG_NEW_CONSTANT(float, 1.41421356237309504880169f, sqrt_2);
      TRNG_NEW_CONSTANT(float, 2.50662827463100050241577f, sqrt_2pi);
      TRNG_NEW_CONSTANT(float, .707106781186547524400845f, one_over_sqrt_2);
      TRNG_NEW_CONSTANT(float, .398942280401432677939946f, one_over_sqrt_2pi)
    };
  
    template<> class constants<double> {
    public:
      TRNG_NEW_CONSTANT(double, 3.14159265358979323846264, pi);
      TRNG_NEW_CONSTANT(double, .915965594177219015054604, catalan);
      TRNG_NEW_CONSTANT(double, 2.71828182845904523536029, e);
      TRNG_NEW_CONSTANT(double, .577215664901532860606512, gamma);
      TRNG_NEW_CONSTANT(double, .693147180559945309417232, ln_2);
      TRNG_NEW_CONSTANT(double, 1.41421356237309504880169, sqrt_2);
      TRNG_NEW_CONSTANT(double, 2.50662827463100050241577, sqrt_2pi);
      TRNG_NEW_CONSTANT(double, .707106781186547524400845, one_over_sqrt_2);
      TRNG_NEW_CONSTANT(double, .398942280401432677939946, one_over_sqrt_2pi)
    };
  
    template<> class constants<long double> {
    public:
      TRNG_NEW_CONSTANT(long double, 3.14159265358979323846264l, pi);
      TRNG_NEW_CONSTANT(long double, .915965594177219015054604l, catalan);
      TRNG_NEW_CONSTANT(long double, 2.71828182845904523536029l, e);
      TRNG_NEW_CONSTANT(long double, .577215664901532860606512l, gamma);
      TRNG_NEW_CONSTANT(long double, .693147180559945309417232l, ln_2);
      TRNG_NEW_CONSTANT(long double, 1.41421356237309504880169l, sqrt_2);
      TRNG_NEW_CONSTANT(long double, 2.50662827463100050241577l, sqrt_2pi);
      TRNG_NEW_CONSTANT(long double, .707106781186547524400845l, one_over_sqrt_2);
      TRNG_NEW_CONSTANT(long double, .398942280401432677939946l, one_over_sqrt_2pi)
    };
    
  }

}

#undef TRNG_NEW_CONSTANT

#endif
