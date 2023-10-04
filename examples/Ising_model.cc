// Copyright (c) 2000-2022, Heiko Bauke
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

#include <iostream>
#include <vector>
#include <queue>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <chrono>
#include <trng/lcg64.hpp>
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
#include <trng/lagfib2xor.hpp>
#include <trng/lagfib4xor.hpp>
#include <trng/uniform_int_dist.hpp>
#include <trng/bernoulli_dist.hpp>

class timer {
private:
  const double _resolution;
  std::chrono::time_point<std::chrono::system_clock> _t;

public:
  void reset() { _t = std::chrono::system_clock::now(); }
  double time() const {
    auto now = std::chrono::system_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(now - _t).count() * 1e-6;
  }
  double resolution() const { return _resolution; }
  timer() : _resolution([]() { return 1e-6; }()), _t(std::chrono::system_clock::now()) {}
};


struct coord {
  int x{0};
  int y{0};
};


class lattice {
private:
  std::vector<int> s;
  int L{0}, L2{0};
  int pos(int) const;

public:
  int size() const;
  void resize(int);
  void fill(int);
  void flipp(coord);
  int get(coord);
  void set(coord, int);
  double energy() const;
  double magnet() const;
  void print() const;
  explicit lattice(int);
  lattice() = default;
};


inline int lattice::pos(int x) const {
  while (x < 0)
    x += L;
  while (x >= L)
    x -= L;
  return x;
}


int lattice::size() const {
  return L;
}


void lattice::resize(int newL) {
  if (newL > 0) {
    L = newL;
    L2 = L * L;
    s.resize(L2);
  } else {
    std::cerr << "negative size\n";
    std::exit(EXIT_FAILURE);
  }
}


inline void lattice::flipp(coord r) {
  s[pos(r.x) + pos(r.y) * L] *= -1;
}


void lattice::fill(int w) {
  for (int i = 0; i < L2; ++i)
    s[i] = w;
}


inline int lattice::get(coord r) {
  return s[pos(r.x) + pos(r.y) * L];
}


void lattice::set(coord r, int w) {
  s[pos(r.x) + pos(r.y) * L] = w;
}


double lattice::energy() const {
  double e{0.0};
  for (int i{0}; i < L; ++i)
    for (int j{0}; j < L; ++j)
      e -= s[i + j * L] * (s[pos(i + 1) + j * L] + s[i + pos(j + 1) * L]);
  return e / L2;
}


double lattice::magnet() const {
  double m{0.0};
  for (int i{0}; i < L; ++i)
    for (int j{0}; j < L; ++j)
      if (s[i + j * L] < 0)
        --m;
      else
        ++m;
  return std::abs(m) / L2;
}


void lattice::print() const {
  for (int i{0}; i < L; ++i) {
    for (int j{0}; j < L; ++j)
      if (s[i + j * L] < 0)
        std::cout << '.';
      else
        std::cout << '#';
    std::cout << '\n';
  }
  std::cout << std::endl;
}


lattice::lattice(int newL) {
  lattice::resize(newL);
}


template<class RNG_type>
void wolffstep(RNG_type &R, lattice &s, double T) {
  std::queue<coord> buffer;
  const double padd = 1.0 - std::exp(-2.0 / T);

  trng::uniform_int_dist U(0, s.size());
  trng::bernoulli_dist<bool> B(padd, true, false);
  coord r;

  r.x = U(R);
  r.y = U(R);
  const int oldspin{s.get(r)};
  s.flipp(r);
  buffer.push(r);
  while (!buffer.empty()) {
    r = buffer.front();
    buffer.pop();
    --r.x;
    if (s.get(r) == oldspin)
      if (B(R)) {
        buffer.push(r);
        s.flipp(r);
      }
    ++r.x;
    ++r.x;
    if (s.get(r) == oldspin)
      if (B(R)) {
        buffer.push(r);
        s.flipp(r);
      }
    --r.x;
    --r.y;
    if (s.get(r) == oldspin)
      if (B(R)) {
        buffer.push(r);
        s.flipp(r);
      }
    ++r.y;
    ++r.y;
    if (s.get(r) == oldspin)
      if (B(R)) {
        buffer.push(r);
        s.flipp(r);
      }
    --r.y;
  }
}


void output(std::vector<double> &Ea, std::vector<double> &ca, int simulations, double E_exact,
            double c_exact);

void output(std::vector<double> &Ea, std::vector<double> &ca, int simulations, double E_exact,
            double c_exact) {
  Ea[simulations] = 0.0;
  ca[simulations] = 0.0;
  Ea[simulations + 1] = 0.0;
  ca[simulations + 1] = 0.0;
  for (int j{0}; j < simulations; ++j) {
    Ea[simulations] += Ea[j] / simulations;
    ca[simulations] += ca[j] / simulations;
  }
  for (int j{0}; j < simulations; ++j) {
    Ea[simulations + 1] += (Ea[j] - Ea[simulations]) * (Ea[j] - Ea[simulations]);
    ca[simulations + 1] += (ca[j] - ca[simulations]) * (ca[j] - ca[simulations]);
  }
  Ea[simulations + 1] /= (simulations - 1.0) * simulations;
  Ea[simulations + 1] = std::sqrt(Ea[simulations + 1]);
  // Student_t(0.99, simulations-1l);
  ca[simulations + 1] /= (simulations - 1.0) * simulations;
  ca[simulations + 1] = std::sqrt(ca[simulations + 1]);
  // Student_t(0.99, simulations-1l);
  std::cout << "\n\t E\t\t c\n";
  for (int j{0}; j < simulations; ++j)
    std::cout << '\t' << std::setprecision(8) << Ea[j] << '\t' << std::setprecision(8) << ca[j]
              << '\n';
  std::cout << "\t--------------\t--------------" << '\n'
            << "mean\t" << std::setprecision(8) << Ea[simulations] << '\t'
            << std::setprecision(8) << ca[simulations] << '\n'
            << "Del\t" << std::setprecision(8) << Ea[simulations] - E_exact << '\t'
            << std::setprecision(8) << ca[simulations] - c_exact << '\n'
            << "sig\t" << std::setprecision(8) << Ea[simulations + 1] << '\t'
            << std::setprecision(8) << ca[simulations + 1] << '\n'
            << "Del/sig\t" << std::setprecision(8)
            << std::abs(Ea[simulations] - E_exact) / Ea[simulations + 1] << '\t'
            << std::setprecision(8) << std::abs(ca[simulations] - c_exact) / ca[simulations + 1]
            << std::endl;
}


template<class RNG_type>
void wolff_main(RNG_type &R, long runs, long split, long L) {
  double E_exact, c_exact;
  switch (L) {
    case 8:
      E_exact = -1.4915891074397066;
      c_exact = 1.1455592398944086;
      break;
    case 12:
      E_exact = -1.4659608164862789;
      c_exact = 1.3529506829072697;
      break;
    case 16:
      E_exact = -1.4530648528134771;
      c_exact = 1.4987049594000261;
      break;
    case 20:
      E_exact = -1.4453094678058525;
      c_exact = 1.6111614949041113;
      break;
    case 24:
      E_exact = -1.4401334960573388;
      c_exact = 1.7027336877232671;
      break;
    case 28:
      E_exact = -1.4364340850836483;
      c_exact = 1.7799744882644384;
      break;
    case 32:
      E_exact = -1.4336584661462483;
      c_exact = 1.8467675900395589;
      break;
    case 36:
      E_exact = -1.4314991053179871;
      c_exact = 1.9056050418011437;
      break;
    case 40:
      E_exact = -1.4297713123073425;
      c_exact = 1.9581816502509387;
      break;
    case 44:
      E_exact = -1.4283574829971357;
      c_exact = 2.0057024700476327;
      break;
    case 48:
      E_exact = -1.4271791793855239;
      c_exact = 2.0490550151069595;
      break;
    case 52:
      E_exact = -1.4261820801536625;
      c_exact = 2.0889118621693695;
      break;
    case 56:
      E_exact = -1.4253273745116323;
      c_exact = 2.1257948956620735;
      break;
    case 60:
      E_exact = -1.4245865955789106;
      c_exact = 2.1601172105639529;
      break;
    case 64:
      E_exact = -1.4239383898330109;
      c_exact = 2.1922113931405711;
      break;
    default:
      std::cerr << "invalid lattice size, try 8, 12, ..., 64\n";
      std::exit(EXIT_FAILURE);
  }

  // const double
  const int simulations{10};
  std::vector<double> Ea(simulations + 2);
  std::vector<double> ca(simulations + 2);

  std::cout.setf(std::ios::fixed);
  std::cout.setf(std::ios::showpoint);
  const double T{2.0 / std::log1p(std::sqrt(2.0))};
  lattice s(L);
  s.fill(-1);
  std::cout << "Generator : " << R.name() << '\n'
            << "Splitting level : " << split << '\n'
            << '\n'
            << "T = " << T << '\n'
            << "Lattice = " << L << "x" << L << '\n'
            << "Samples = " << runs << '\n';

  for (int i{0}; i < 2 * runs; ++i)
    wolffstep(R, s, T);
  timer t;
  for (int j{0}; j < simulations; ++j) {
    double E{0.0};
    double E2{0.0};
    for (int i{0}; i < runs; ++i) {
      wolffstep(R, s, T);
      const double q{s.energy()};
      E += q;
      E2 += q * q;
    }
    E /= runs;
    E2 /= runs;
    double c = (L * L) / (T * T) * (E2 - E * E);
    Ea[j] = E;
    ca[j] = c;
  }
  output(Ea, ca, simulations, E_exact, c_exact);
  std::cout << '\n' << "Time: " << t.time() << " sec." << std::endl;
}

int main(int argc, char *argv[]) {
  long runs{0}, split{0}, L{0};
  std::string generator;
  int argi{1};
  if (argc == 1) {
    std::cerr << "Tina's Random Number Generator Library\n\n"
              << "(P) & (C) by Heiko Bauke, 2000-2022\n\n"
              << "two-domensional Ising model (Wolff algorithm)\n"
              << "---------------------------------------------\n\n"
              << "synopsis:\n"
              << "$ " << argv[0] << " --gen generator --runs runs --split split --size size\n"
              << "try:\n"
              << "$ " << argv[0] << " --gen lcg64 --runs 100000 --split 1 --size 16\n";
    std::exit(EXIT_FAILURE);
  }
  while (argi < argc) {
    std::string arg(argv[argi]);
    if (arg == "--gen") {
      ++argi;
      if (argi < argc) {
        generator = argv[argi];
        ++argi;
      } else {
        std::cerr << "missung argument for parameter " << arg << '\n';
        std::exit(EXIT_FAILURE);
      }
    } else if (arg == "--runs") {
      ++argi;
      if (argi < argc) {
        runs = std::atoi(argv[argi]);
        ++argi;
      } else {
        std::cerr << "missung argument for parameter " << arg << '\n';
        std::exit(EXIT_FAILURE);
      }
    } else if (arg == "--split") {
      ++argi;
      if (argi < argc) {
        split = std::atoi(argv[argi]);
        ++argi;
      } else {
        std::cerr << "missung argument for parameter " << arg << '\n';
        std::exit(EXIT_FAILURE);
      }
    } else if (arg == "--size") {
      ++argi;
      if (argi < argc) {
        L = std::atoi(argv[argi]);
        ++argi;
      } else {
        std::cerr << "missung argument for parameter " << arg << '\n';
        std::exit(EXIT_FAILURE);
      }
    } else {
      std::cerr << "unknown argument " << arg << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
  {
    trng::lcg64 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg2 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg3 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg3s r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg4 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg5 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::mrg5s r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn2 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn3 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn3s r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn4 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn5 r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::yarn5s r;
    r.split(split, 0);
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::r250_32 r;
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  {
    trng::Ziff_32 r;
    if (generator == r.name())
      wolff_main(r, runs, split, L);
  }
  return EXIT_SUCCESS;
}
