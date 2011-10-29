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

#include <iostream>
#include <vector>
#include <queue>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
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
#include <trng/uniform_int_dist.hpp>
#include <trng/bernoulli_dist.hpp>

#if defined __unix__
# include <unistd.h>
# include <sys/time.h>
# include <sys/times.h>
#else
# include <ctime>
#endif

class timer {
private:
  const double _resolution;
  double _t;
  double get_time() {
#if defined __unix__
    struct timeval  tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return static_cast<double>(tv.tv_sec)+static_cast<double>(tv.tv_usec)*1e-6;
#else
    return static_cast<double>(std::clock())*_resolution;
#endif
  }
public:
  void reset() { _t=get_time(); }
  double time() { return get_time()-_t; }
  double resolution() const { return _resolution; }
  timer() :
#if defined __unix__
    _resolution(1e-6),
#else
    _resolution(1.0/CLOCKS_PER_SEC),
#endif
    _t(get_time()) { }
};


typedef struct {
  int x;
  int y;
} koord;

class lattice {
private:
  std::vector<int> s;
  int L, L2;
  int pos(int);
public:
  int size(void);
  void resize(int, int);
  void fill(int);
  void flipp(koord);
  int get(koord);
  void set(koord, int);
  double energy(void);
  double magnet(void);
  void print(void);
  lattice(int);
  lattice();
};

inline int lattice::pos(int x) {
  while (x<0)
    x+=L;
  while (x>=L)
    x-=L;
  return x;
}

int lattice::size(void) {
  return L;
}

void lattice::resize(int newL, int w=-1) {
  if (newL>0) {
    L=newL;
    L2=L*L;
    s.resize(L2);
  } else {
    std::cerr << "negative size\n";
    std::exit(EXIT_FAILURE);
  }
}

inline void lattice::flipp(koord r) {
   s[pos(r.x)+pos(r.y)*L]*=-1;;
}

void lattice::fill(int w) {
  for (int i=0; i<L2; ++i)
    s[i]=w;
}

inline int lattice::get(koord r) {
  return s[pos(r.x)+pos(r.y)*L];
}

void lattice::set(koord r, int w) {
  s[pos(r.x)+pos(r.y)*L]=w;
}

double lattice::energy(void) {
  double e=0.0;
  for (int i=0; i<L; ++i)
    for (int j=0; j<L; ++j)
      e-=s[i+j*L]*(s[pos(i+1)+j*L]+s[i+pos(j+1)*L]);
  return e/L2;
}

double lattice::magnet(void) {
  double m=0.0;
  for (int i=0; i<L; ++i)
    for (int j=0; j<L; ++j)
      if (s[i+j*L]<0)
	--m;
      else
	++m;
  return std::abs(m)/L2;
}

void lattice::print(void) {
  int i, j;
  for (int i=0; i<L; ++i) {
    for (int j=0; j<L; ++j)
      if (s[i+j*L]<0)
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

lattice::lattice() {
  L=0;
  L2=0;
}


template<class RNG_type>
void wolffstep(RNG_type &R, lattice &s, double T) {
  std::queue<koord> buffer;
  const double padd=1.0-std::exp(-2.0/T);;
  int oldspin;
  trng::uniform_int_dist U(0, s.size());
  trng::bernoulli_dist<bool> B(padd, true, false);
  koord r;

  r.x=U(R);
  r.y=U(R);
  oldspin=s.get(r);
  s.flipp(r);
  buffer.push(r);
  while (!buffer.empty()) {
    r=buffer.front();
    buffer.pop();
    --r.x;
    if (s.get(r)==oldspin)
      if (B(R)) {
	buffer.push(r);
	s.flipp(r);
      }
    ++r.x;
    ++r.x;
    if (s.get(r)==oldspin)
      if (B(R)) {
	buffer.push(r);
	s.flipp(r);
      }
    --r.x;
    --r.y;
    if (s.get(r)==oldspin)
      if (B(R)) {
	buffer.push(r);
	s.flipp(r);
      }
    ++r.y;
    ++r.y;
    if (s.get(r)==oldspin)
      if (B(R)) {
	buffer.push(r);
	s.flipp(r);
      }
    --r.y;
  }
}

void output(std::vector<double> &Ea, std::vector<double> &ca, 
	    int simulations, double E_exact, double c_exact) {
  int j;
  Ea[simulations]=0.0;
  ca[simulations]=0.0;
  Ea[simulations+1]=0.0;
  ca[simulations+1]=0.0;
  for (j=0; j<simulations; ++j) {
    Ea[simulations]+=Ea[j]/simulations;
    ca[simulations]+=ca[j]/simulations;
  }
  for (j=0; j<simulations; ++j) {
    Ea[simulations+1]+=(Ea[j]-Ea[simulations])*(Ea[j]-Ea[simulations]);
    ca[simulations+1]+=(ca[j]-ca[simulations])*(ca[j]-ca[simulations]);
  }
  Ea[simulations+1]/=(simulations-1.0)*simulations;
  Ea[simulations+1]=std::sqrt(Ea[simulations+1]);
  //Student_t(0.99, simulations-1l);
  ca[simulations+1]/=(simulations-1.0)*simulations;
  ca[simulations+1]=std::sqrt(ca[simulations+1]);
  //Student_t(0.99, simulations-1l);
  std::cout << "\n\t E\t\t c\n";
  for (j=0; j<simulations; ++j)
    std::cout << '\t' << std::setprecision(8) << Ea[j] << '\t' 
	      << std::setprecision(8) << ca[j] << '\n';
  std::cout << "\t--------------\t--------------" << '\n'
	    << "mean\t" << std::setprecision(8) << Ea[simulations] << '\t' 
	    << std::setprecision(8) << ca[simulations] << '\n'
	    << "Del\t" 
	    << std::setprecision(8) << Ea[simulations]-E_exact << '\t' 
	    << std::setprecision(8) << ca[simulations]-c_exact << '\n'
	    << "sig\t" << std::setprecision(8) << Ea[simulations+1] << '\t' 
	    << std::setprecision(8) << ca[simulations+1] << '\n'
	    << "Del/sig\t" 
	    << std::setprecision(8) 
	    << std::abs(Ea[simulations]-E_exact)/Ea[simulations+1] << '\t' 
	    << std::setprecision(8) 
	    << std::abs(ca[simulations]-c_exact)/ca[simulations+1] << std::endl;
}


template<class RNG_type>
void wolff_main(RNG_type &R, long runs, long split) {
  const int L=16;
  const double E_exact=-1.453064854; // 16x16
  const double c_exact= 1.498704952; // 16x16
  // const double E_exact=-1.465960817; // 12x12
  // const double c_exact= 1.352950680; // 12x12
  // const double E_exact=-1.491589107;  // 8x8
  // const double c_exact= 1.145559238;  // 8x8
  const int simulations=10;
  int i, j;
  double T, E, E2, q, c;
  std::vector<double> Ea(simulations+2);
  std::vector<double> ca(simulations+2);

  std::cout.setf(std::ios::fixed);
  std::cout.setf(std::ios::showpoint);
  T=2.0/std::log(1.0+std::sqrt(2.0));
  lattice s(L);
  s.fill(-1);
  R.split(split, 0);
  std::cout << "Generator : " << R.name() << '\n'
	    << "Splitting level : " << split << '\n'
	    << '\n'
	    << "T = " << T << '\n'
	    << "Lattice = " << L << "x" << L << '\n'
	    << "Samples = " << runs << '\n';

  for (i=0; i<2*runs; ++i) 
    wolffstep(R, s, T);
  timer t;  
  for (j=0; j<simulations; ++j) { 
    E=0.0;
    E2=0.0;
    for (i=0; i<runs; ++i) {
      wolffstep(R, s, T);
      q=s.energy();
      E+=q;
      E2+=q*q;
    }
    E/=runs;
    E2/=runs;
    c=(L*L)/(T*T)*(E2-E*E);
    Ea[j]=E;
    ca[j]=c;
  }
  output(Ea, ca, simulations, E_exact, c_exact);
  std::cout << '\n' 
	    << "Time: " << t.time() << " sec." << std::endl;
}

int main (int argc, char *argv[]) {
  long runs, split;
  std::string generator;
  if (argc!=4) {
    std::cerr << "Tina's Random Number Generator Library\n\n"
	      << "(P) & (C) by Heiko Bauke, Magdeburg 2001-2006\n\n"
	      << "two-domensional Ising model (Wolff algorithm)\n"
	      << "---------------------------------------------\n\n"
	      << "synopsis:\n"
	      << "$ " << argv[0] << " generator runs split\n"
	      << "usefull arguments: generator=lcg64; runs=100000; split=1\n";
    exit(EXIT_FAILURE);
  }
  generator=argv[1];
  runs=atoi(argv[2]);
  split=atoi(argv[3]);
  { trng::lcg64 r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg2 r;   if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg3 r;   if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg3s r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg4 r;   if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg5 r;   if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::mrg5s r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn2 r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn3 r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn3s r; if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn4 r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn5 r;  if (generator==r.name()) wolff_main(r, runs, split); }
  { trng::yarn5s r; if (generator==r.name()) wolff_main(r, runs, split); }
  return EXIT_SUCCESS;
}
