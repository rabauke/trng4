// Copyright (C) 2000-2010 Heiko Bauke <heiko.bauke@mpi-hd.mpg.de>
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

#include <cassert>
#include <cstdio>
#include <istream>
#include <ostream>
#include <iomanip>
#include <ios>
#include <cstring>
#include <vector>
#include <iterator>
#include <trng/limits.hpp>
#include <trng/uniformxx.hpp>

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
	while (i<len and !(in.get(c) and c!=d.str[i])) {
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
	  if (c==EOF or !(c==' ' or c=='\t' or c=='\n'))
	    break;
	  in.get();
	}
	return in;
      }
    }; 
    
    inline ignore_spaces_cl ignore_spaces() {
      return ignore_spaces_cl();
    }
    
    // -------------------------------------------------------------------
      
    template<typename T>
    class io_range {
    private:
      T first, last;
      const char *delim_str;
    public:
      io_range(const T first, const T last,
	       const char *delim_str=0) : first(first),
					  last(last),
					  delim_str(delim_str) {
      }

      template<typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, const io_range &IO_range) {
	T pos(IO_range.first);
	while (out and pos!=IO_range.last) {
	  out << (*pos);
	  ++pos;
	  if (pos!=IO_range.last and IO_range.delim_str!=0)
	    out << IO_range.delim_str;
	}
	return out;
      }
    
      template<typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> &in, const io_range &IO_range) {
	T pos(IO_range.first);
	while (in and pos!=IO_range.last) {
	  in >> (*pos);
	  ++pos;
	  if (pos!=IO_range.last)
	    in >> delim(IO_range.delim_str);
	}
	return in;
      }
    
    };
    
    template<typename T>
    inline io_range<T> make_io_range(T first, T last, const char *delim_str=0) {
      return io_range<T>(first, last, delim_str);
    }

    // ---------------------------------------------------------------

    template<long m>
    struct log2_floor {
      enum { value = 1 + log2_floor<m/2>::value };
    };

    template<>
    struct log2_floor<0> { enum { value = 0 }; };

    template<>
    struct log2_floor<1> { enum { value = 0 }; };
    
    template<long m>
    struct log2_ceil {
      enum { value = (1ul<<log2_floor<m>::value) < m ? log2_floor<m>::value+1 : log2_floor<m>::value }; 
    };

    // ---------------------------------------------------------------

    template<long m, long r> 
    class modulo_helper;

    template<long m> 
    class modulo_helper<m, 0> {
      static const long e=log2_ceil<m>::value;
      static const long k=(1ul<<e)-m;
      static const unsigned long mask=(1ul<<e)-1ul;
    public:
      inline static long modulo(unsigned long long x) {
	if (mask==m) {
	  unsigned long y=(x&mask)+(x>>e);
	  if (y>=m)
	    y-=m;
	  return y;
	} else if (static_cast<long long>(k)*static_cast<long long>(k+2)<=m) {
	  x=(x&mask)+(x>>e)*k;
	  x=(x&mask)+(x>>e)*k;
	  if (x>=m)
	    x-=m;
	  return x;
	} else {
	  return x%m;
	}
      }
    };

    template<long m> 
    class modulo_helper<m, 1> {
      static const long e=log2_ceil<m>::value;
      static const long k=(1ul<<e)-m;
      static const unsigned long mask=(1ul<<e)-1ul;
    public:
      inline static long modulo(unsigned long long x) {
	if (mask==m) {
	  unsigned long long y=(x&mask)+(x>>e);
	  if (y>=2ull*m)
	    y-=2ull*m;
	  if (y>=m)
	    y-=m;
	  return y;
	} else if (static_cast<long long>(k)*static_cast<long long>(k+2)<=m) {
	  x=(x&mask)+(x>>e)*k;
	  x=(x&mask)+(x>>e)*k;
	  if (x>=2ull*m) x-=2ull*m; 
	  if (x>=m)
	    x-=m;
	  return x;
	} else {
	  return x%m;
	}
      }
    };

    template<long m> 
    class modulo_helper<m, 2> {
      static const long e=log2_ceil<m>::value;
      static const long k=(1ul<<e)-m;
      static const unsigned long mask=(1ul<<e)-1ul;
    public:
      inline static long modulo(unsigned long long x) {
	if (mask==m) {
	  unsigned long long y=(x&mask)+(x>>e);
	  if (y>=4ull*m) y-=4ull*m;
	  if (y>=2ull*m) y-=2ull*m;
	  if (y>=m)
	    y-=m;
	  return y;
	} else if (static_cast<long long>(k)*static_cast<long long>(k+2)<=m) {
	  x=(x&mask)+(x>>e)*k;
	  x=(x&mask)+(x>>e)*k;
	  if (x>=4ull*m) x-=4ull*m; 
	  if (x>=2ull*m) x-=2ull*m; 
	  if (x>=m)
	    x-=m;
	  return x;
	} else {
	  return x%m;
	}
      }
    };
    
    template<long m, long r> 
    inline long modulo(unsigned long long x) {
      return modulo_helper<m, log2_floor<r>::value >::modulo(x);
    }
    
    // ---------------------------------------------------------------

    template<long m, long b>
    class power {
      unsigned long b_power0[0x10000l], b_power1[0x08000l];

      inline long pow(long n) {
        long long p(1ll), t(b);
        while (n>0) {
          if ((n&0x1)==0x1)
            p=(p*t)%m;
          t=(t*t)%m;
          n/=2;
        }
        return static_cast<long>(p);
      }
      
      power & operator=(const power &);
      power(const power &);
      
    public:
      power() {
        for (long i(0l); i<0x10000l; ++i)
          b_power0[i]=pow(i);
        for (long i(0l); i<0x08000l; ++i)
          b_power1[i]=pow(i*0x10000l);
      }
      inline long operator()(long n) const {
        return modulo<m, 1>(static_cast<unsigned long long>(b_power1[n>>16])*
			    static_cast<unsigned long long>(b_power0[n&0xffff]));
      }

    };

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
        
  }  
  
}

#endif
