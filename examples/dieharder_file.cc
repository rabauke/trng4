// Copyright (c) 2000-2026, Heiko Bauke
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
#include <string>
#include <sstream>
#include <functional>
#include <map>
#if defined _MSC_VER && __cplusplus <= 201703
#include <ciso646>
#endif
#include <trng/uniform01_dist.hpp>
#include <trng/lcg64.hpp>
#include <trng/lcg64_shift.hpp>
#include <trng/lcg64_count_shift.hpp>
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
#include "trng/mt19937.hpp"
#include "trng/mt19937_64.hpp"
#include "trng/count128_lcg_shift.hpp"


template<typename T>
bool try_parse(const char *const str, T &res) {
  std::stringstream stream;
  stream << str;
  res = T{};
  stream >> res;
  return stream.eof() and not stream.fail();
}


template<typename R>
void generate(std::size_t samples, unsigned long seed) {
  R r;
  std::cout << "#==================================================================\n"
            << "# generator " << r.name() << "  seed = " << seed << "\n"
            << "#==================================================================\n"
            << "type: d\n"
            << "count: " << samples << "\n"
            << "numbit: " << trng::int_math::log2_ceil(r.max() - r.min()) << '\n';
  for (std::size_t j{0}; j < samples; ++j)
    std::cout << r() - r.min() << '\n';
}


template<typename R>
void add_generator(
    std::map<std::string, std::function<void(std::size_t, unsigned long)>> &func_map) {
  const R r;
  func_map[r.name()] = [&](std::size_t samples, unsigned long seed) {
    generate<R>(samples, seed);
  };
}


class command_line_argument_error : public std::runtime_error {
  using base = std::runtime_error;

public:
  using base::base;
};


int main(int argc, char *argv[]) {
  std::map<std::string, std::function<void(std::size_t, unsigned long)>> func_map;

  add_generator<trng::lcg64>(func_map);
  add_generator<trng::lcg64_shift>(func_map);
  add_generator<trng::lcg64_count_shift>(func_map);
  add_generator<trng::mrg2>(func_map);
  add_generator<trng::mrg3>(func_map);
  add_generator<trng::mrg3s>(func_map);
  add_generator<trng::mrg4>(func_map);
  add_generator<trng::mrg5>(func_map);
  add_generator<trng::mrg5s>(func_map);
  add_generator<trng::yarn2>(func_map);
  add_generator<trng::yarn3>(func_map);
  add_generator<trng::yarn3s>(func_map);
  add_generator<trng::yarn4>(func_map);
  add_generator<trng::yarn5>(func_map);
  add_generator<trng::yarn5s>(func_map);
  add_generator<trng::mt19937>(func_map);
  add_generator<trng::mt19937_64>(func_map);
  add_generator<trng::count128_lcg_shift>(func_map);

  try {
    if (argc != 3 and argc != 4)
      throw command_line_argument_error("wrong number of arguments");
    const auto generate_iter{func_map.find(argv[1])};
    if (generate_iter == func_map.end())
      throw command_line_argument_error("unknown generator");
    std::size_t samples{0};
    if (not try_parse(argv[2], samples))
      throw command_line_argument_error("illegal number of samples");
    unsigned long seed{0};
    if (argc == 4 and not try_parse(argv[3], seed))
      throw command_line_argument_error("illegal seed value");

    (generate_iter->second)(samples, seed);
  } catch (const command_line_argument_error &ex) {
    std::cerr
        << "error: " << ex.what() << "\n\n"
        << "Create an input file for the Dieharder Random Number Test Suite that can be\n"
        << "used as input for generator \"202\", see Dieharder documentation for details.\n\n"
        << "- https://webhome.phy.duke.edu/~rgb/General/dieharder.php\n"
        << "- https://github.com/eddelbuettel/dieharder\n\n"
        << "Usage :\n"
        << argv[0] << " [PRNG] [number of samples] [seed]\n\n"
        << "List of possible PRNGs:\n";
    for (auto i{func_map.begin()}; i != func_map.end(); ++i)
      std::cerr << "  " << i->first << '\n';
  }
  return EXIT_SUCCESS;
}
