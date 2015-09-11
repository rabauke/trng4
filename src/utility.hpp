// Copyright (c) 2000-2015, Heiko Bauke
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
#include <trng/cuda.hpp>
#include <trng/limits.hpp>
#include <trng/uniformxx.hpp>

namespace trng {
  
  namespace utility {
    
    template<typename T>
    TRNG_CUDA_ENABLE
    inline const T & min(const T &a, const T &b) {
      return a<=b ? a : b;
    }

    template<typename T>
    TRNG_CUDA_ENABLE
    inline const T & max(const T &a, const T &b) {
      return a>=b ? a : b;
    }
    
    // ---------------------------------------------------------------

    template<typename T>
    TRNG_CUDA_ENABLE
    inline void swap(T &a, T &b) {
      T c(b);
      b=a;
      a=c;
    }

    // ---------------------------------------------------------------

    template<typename T>
    inline void throw_this(const T &x) {
      throw x;
    }

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
    class io_range;

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

      template<typename S, typename char_t, typename traits_t>
      friend std::basic_ostream<char_t, traits_t> & 
      operator<<(std::basic_ostream<char_t, traits_t> &out, const io_range<S> &IO_range);
    
      template<typename S, typename char_t, typename traits_t>
      friend std::basic_istream<char_t, traits_t> & 
      operator>>(std::basic_istream<char_t, traits_t> &in, const io_range<S> &IO_range);
    
    };

    template<typename T, typename char_t, typename traits_t>
    std::basic_ostream<char_t, traits_t> & 
    operator<<(std::basic_ostream<char_t, traits_t> &out, const io_range<T> &IO_range) {
      T pos(IO_range.first);
      while (out and pos!=IO_range.last) {
	out << (*pos);
	++pos;
	if (pos!=IO_range.last and IO_range.delim_str!=0)
	    out << IO_range.delim_str;
      }
      return out;
    }
    
    template<typename T, typename char_t, typename traits_t>
    std::basic_istream<char_t, traits_t> & 
    operator>>(std::basic_istream<char_t, traits_t> &in, const io_range<T> &IO_range) {
      T pos(IO_range.first);
      while (in and pos!=IO_range.last) {
	in >> (*pos);
	++pos;
	if (pos!=IO_range.last)
	  in >> delim(IO_range.delim_str);
      }
      return in;
    }

    template<typename T>
    inline io_range<T> make_io_range(T first, T last, const char *delim_str=0) {
      return io_range<T>(first, last, delim_str);
    }
    
    // ---------------------------------------------------------------
    
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
