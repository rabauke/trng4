// Copyright (C) 2006 Heiko Bauke <heiko.bauke@physik.uni-magdeburg.de>
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

#if !(defined TRNG_UTILITY_HPP)

#define TRNG_UTILITY_HPP

#include <istream>
#include <ostream>
#include <iomanip>
#include <ios>
#include <cstring>
#include <vector>
#include <iterator>
#include <trng/limits.hpp>

namespace trng {
  
  namespace utility {
    
    class delim_str;
    class delim_c;
    inline delim_str delim(const char * const);
    inline delim_c delim(char);
    
    class ignore_spaces_cl;
    inline ignore_spaces_cl ignore_spaces();

    long modulo_invers(long a, long m);
    void gauss(std::vector<long> &a,
	       std::vector<long> &b, long m);
    void matrix_mult(const std::vector<long> &a,
		     const std::vector<long> &b,
		     std::vector<long> &c, long m);
    void matrix_vec_mult(const std::vector<long> &a,
			 const std::vector<long> &b,
			 std::vector<long> &c, long m);
    
    template<typename iter>
    inline int discrete(double x, iter, iter);
    
    // ---------------------------------------------------------------

    class delim_str {
      const char * str;
    public:
      delim_str(const char * str_) : str(str_) { }
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> & in, 
		 const delim_str &d) {
	char c;
	std::size_t len(std::strlen(d.str)), i(0);
	while (i<len && !(in.get(c) && c!=d.str[i])) {
	  ++i;
	}
 	if (i<len)
 	  in.setstate(std::ios::failbit);
	return in;
      }
    };
    
    class delim_c {
      const char c;
    public:
      delim_c(char c_) : c(c_) { }
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> & in, 
		 const delim_c &d) {
 	char c;
	in.get(c);
 	if (c!=d.c)
 	  in.setstate(std::ios::failbit);
       	return in;
      }
    };
    
    inline delim_str delim(const char * const str) {
      return delim_str(str);
    }
    
    inline delim_c delim(char c) {
      return delim_c(c);
    }  

    // -----------------------------------------------------------------

    class ignore_spaces_cl {
    public:
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> & in, 
		 const ignore_spaces_cl &) {
	while (true) {
	  int c(in.peek());
	  if (c==EOF || !(c==' ' || c=='\t' || c=='\n'))
	    break;
	  in.get();
	}
	return in;
      }
    }; 
    
    inline ignore_spaces_cl ignore_spaces() {
      return ignore_spaces_cl();
    }
    
    // ---------------------------------------------------------------

    template<long m> 
    inline long modulo(unsigned long long x) {
      return x%m;
    }
    
    // 2^31 - 1 
    template<>
    inline long modulo<2147483647l>(unsigned long long x) {
      unsigned long y(static_cast<unsigned long>(x&0x7ffffffful)+
		      static_cast<unsigned long>(x>>31));
      return (y>=2147483647ul) ? (y-2147483647ul) : y;
    }

    // 2^31 - 21069
    template<>
    inline long modulo<2147462579l>(unsigned long long x) {
      x=(x&0x7ffffffful)+(x>>31)*21069;
      x=(x&0x7ffffffful)+(x>>31)*21069;
      return (x>=2147462579ul) ? (x-2147462579ul) : x;
    }

    // 2^31 - 22641
    template<>
    inline long modulo<2147461007l>(unsigned long long x) {
      x=(x&0x7ffffffful)+(x>>31)*22641;
      x=(x&0x7ffffffful)+(x>>31)*22641;
      return (x>=2147461007ul) ? (x-2147461007ul) : x;
    }
    
    // ---------------------------------------------------------------

    template<long modulus>
    class power {
      long b;
      unsigned long b_power0[0x10000], b_power1[0x08000];
      
      inline long pow(long n) {
	long long p(1ll), t(b);
	while (n>0) {
	  if ((n&0x1)==0x1)
	    p=(p*t)%modulus;
	t=(t*t)%modulus;
	n/=2;
	}
	return static_cast<long>(p);
      }
      
      void calc_b_power() {
	for (long i(0l); i<0x10000l; ++i)
	  b_power0[i]=pow(i);
	for (long i(0l); i<0x08000l; ++i)
	  b_power1[i]=pow(i*0x10000l);
      }

    public:
      inline long operator()(long n) const {
	return modulo<modulus>(static_cast<unsigned long long>(b_power1[n>>16])*
			       static_cast<unsigned long long>(b_power0[n&0xffff]));
      }
      
      long operator()() const {
	return b;
      }

      power(const long b=2) : b(b) {
	calc_b_power();
      }
      
      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, const power &g) {
	out << g.b;
	return out;
      }

      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> &in, power &g) {
	long b, modulus;
	in >> b;
	if (in) 
	  g=power(b);
	return in;
      }
    };
    
    template<long modulus>
    inline bool operator==(const power<modulus> &p1, const power<modulus> &p2) {
      return p1()==p2();
    }

    template<long modulus>
    inline bool operator!=(const power<modulus> &p1, const power<modulus> &p2) {
      return p1()!=p2();
    }
    
    // -----------------------------------------------------------------

    template <unsigned int x> 
    struct ceil2 {
      enum { result=2*ceil2<(x+1)/2>::result };
    };
    
    template <>
    struct ceil2<0> {
      enum { result=0 };
    };
    
    template <>
    struct ceil2<1> {
      enum { result=1 };
    };
    
    template <unsigned int x> 
    struct mask {
      enum { result=ceil2<x+1>::result-1 };
    };

    // -----------------------------------------------------------------

    // random number
    template<typename iter>
    inline int discrete(double x, iter first, iter last) {
      typedef typename std::iterator_traits<iter>::difference_type
	difference_type;
      if (x<(*first))
	return 0;
      difference_type i1(0), i2(last-first-1);
      while (i2-i1>difference_type(1)) {
	difference_type i3((i2+i1)/2);
	if (x<=first[i3])
	  i2=i3;
	else
	  i1=i3;
      }
      return i2;
    }
    
    // -----------------------------------------------------------------

    template<typename R>
    inline double uniformcc_impl(R &r, long) {
      return static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
    }
    
    template<typename R>
    inline double uniformcc_impl(R &r, unsigned long) {
      return static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
    }

    template<typename R>
    inline double uniformcc_impl(R &r, long long) {
      return static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
    }

    template<typename R>
    inline double uniformcc_impl(R &r, unsigned long long) {
      return static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
    }

    template<typename R>
    inline double uniformcc(R &r) {
      return uniformcc_impl(r, typename R::result_type());
    }

    // ---------------------------------------------------------------------

    template<typename R>
    inline double uniformoo_impl(R &r, long) {
      return (1.0+static_cast<double>(r()-R::min))/
	(2.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformoo_impl(R &r, unsigned long) {
      return (1.0+static_cast<double>(r()-R::min))/
	(2.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformoo_impl(R &r, long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==0.0)
	return math::numeric_limits<double>::epsilon();
      if (t==1.0)
	return 1.0-math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformoo_impl(R &r, unsigned long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==0.0)
	return math::numeric_limits<double>::epsilon();
      if (t==1.0)
	return 1.0-math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformoo(R &r) {
      return uniformoo_impl(r, typename R::result_type());
    }

    // -----------------------------------------------------------------

    template<typename R>
    inline double uniformco_impl(R &r, long) {
      return static_cast<double>(r()-R::min)/
	(1.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformco_impl(R &r, unsigned long) {
      return static_cast<double>(r()-R::min)/
	(1.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformco_impl(R &r, long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==1.0)
	return 1.0-math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformco_impl(R &r, unsigned long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==1.0)
	return 1.0-math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformco(R &r) {
      return uniformco_impl(r, typename R::result_type());
    }

    // -----------------------------------------------------------------

    template<typename R>
    inline double uniformoc_impl(R &r, long) {
      return (1.0+static_cast<double>(r()-R::min))/
	(1.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformoc_impl(R &r, unsigned long) {
      return static_cast<double>(r()-R::min)/
	(1.0+static_cast<double>(R::max-R::min));
    }

    template<typename R>
    inline double uniformoc_impl(R &r, long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==0.0)
	return math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformoc_impl(R &r, unsigned long long) {
      double t=static_cast<double>(r()-R::min)/
	static_cast<double>(R::max-R::min);
      if (t==0.0)
	return math::numeric_limits<double>::epsilon();
      return t;
    }

    template<typename R>
    inline double uniformoc(R &r) {
      return uniformoc_impl(r, typename R::result_type());
    }
    
  }
  
}

#endif
