// Copyright (c) 2000-2024, Heiko Bauke
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

#include <cstdlib>
#include <iostream>
#include <exception>
#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/mrg2.hpp>
#include <trng/mrg3.hpp>
#include <trng/mrg3s.hpp>
#include <trng/mrg4.hpp>
#include <trng/mrg5.hpp>
#include <trng/mrg5s.hpp>
#include <trng/yarn2.hpp>
#include <trng/yarn3.hpp>
#include <trng/yarn3s.hpp>
#include <trng/yarn4.hpp>
#include <trng/yarn5.hpp>
#include <trng/yarn5s.hpp>
#include <trng/mt19937.hpp>
#include <trng/mt19937_64.hpp>
#include <trng/lagfib2xor.hpp>
#include <trng/lagfib2plus.hpp>
#include <trng/lagfib4xor.hpp>
#include <trng/lagfib4plus.hpp>


template<typename R>
void sample_output(R &r, std::string name) {
  while (name.length() < 32)
    name += ' ';
  std::cout << name;
  for (int i{0}; i < 15; ++i)
    std::cout << r() << '\t';
  std::cout << r() << '\n';
}


int main() {
  try {
    {
      trng::lcg64 r;
      sample_output(r, "trng::lcg64");
    }
    {
      trng::lcg64_shift r;
      sample_output(r, "trng::lcg64_shift");
    }
    {
      trng::mrg2 r;
      sample_output(r, "trng::mrg2");
    }
    {
      trng::mrg3 r;
      sample_output(r, "trng::mrg3");
    }
    {
      trng::mrg3s r;
      sample_output(r, "trng::mrg3s");
    }
    {
      trng::mrg4 r;
      sample_output(r, "trng::mrg4");
    }
    {
      trng::mrg5 r;
      sample_output(r, "trng::mrg5");
    }
    {
      trng::mrg5s r;
      sample_output(r, "trng::mrg5s");
    }
    {
      trng::yarn2 r;
      sample_output(r, "trng::yarn2");
    }
    {
      trng::yarn3 r;
      sample_output(r, "trng::yarn3");
    }
    {
      trng::yarn3s r;
      sample_output(r, "trng::yarn3s");
    }
    {
      trng::yarn4 r;
      sample_output(r, "trng::yarn4");
    }
    {
      trng::yarn5 r;
      sample_output(r, "trng::yarn5");
    }
    {
      trng::yarn5s r;
      sample_output(r, "trng::yarn5s");
    }
    {
      trng::mt19937 r;
      sample_output(r, "trng::mt19937");
    }
    {
      trng::mt19937_64 r;
      sample_output(r, "trng::mt19937_64");
    }
    {
      trng::lagfib2xor_19937_64 r;
      sample_output(r, "trng::lagfib2xor_19937_64");
    }
    {
      trng::lagfib4xor_19937_64 r;
      sample_output(r, "trng::lagfib4xor_19937_64");
    }
    {
      trng::lagfib4plus_19937_64 r;
      sample_output(r, "trng::lagfib2plus_19937_64");
    }
    {
      trng::lagfib2plus_19937_64 r;
      sample_output(r, "trng::lagfib4plus_19937_64");
    }
  } catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
  } catch (...) {
    std::cerr << "something went wrong" << std::endl;
  }
  return EXIT_SUCCESS;
}
