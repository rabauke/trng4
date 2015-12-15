AUTOMAKE_OPTIONS = gnu

lib_LTLIBRARIES = libtrng4.la

AM_CPPFLAGS = -I$(top_builddir)
AM_CXXFLAGS = @BOOST_CXXFLAGS@

libtrng4_la_SOURCES = minstd.cc \
		      lcg64.cc lcg64_shift.cc \
		      mt19937.cc mt19937_64.cc \
		      mrg2.cc mrg3.cc mrg3s.cc mrg4.cc mrg5.cc mrg5s.cc \
		      yarn2.cc yarn3.cc yarn3s.cc yarn4.cc yarn5.cc yarn5s.cc

pkginclude_HEADERS = config.hpp cuda.hpp utility.hpp uniformxx.hpp static_assertion.hpp \
		     generate_canonical.hpp minstd.hpp \
		     limits.hpp math.hpp int_types.hpp int_math.hpp special_functions.hpp constants.hpp \
		     lcg64.hpp lcg64_shift.hpp \
		     mt19937.hpp mt19937_64.hpp \
		     mrg2.hpp mrg3.hpp mrg3s.hpp mrg4.hpp mrg5.hpp mrg5s.hpp \
		     yarn2.hpp yarn3.hpp yarn3s.hpp yarn4.hpp yarn5.hpp yarn5s.hpp \
		     lagfib2xor.hpp lagfib4xor.hpp lagfib2plus.hpp lagfib4plus.hpp \
		     uniform_dist.hpp uniform01_dist.hpp uniform_int_dist.hpp \
		     exponential_dist.hpp twosided_exponential_dist.hpp \
		     normal_dist.hpp truncated_normal_dist.hpp correlated_normal_dist.hpp \
                     maxwell_dist.hpp \
		     cauchy_dist.hpp logistic_dist.hpp lognormal_dist.hpp pareto_dist.hpp \
		     powerlaw_dist.hpp tent_dist.hpp weibull_dist.hpp \
		     extreme_value_dist.hpp gamma_dist.hpp beta_dist.hpp chi_square_dist.hpp student_t_dist.hpp snedecor_f_dist.hpp\
		     rayleigh_dist.hpp \
		     bernoulli_dist.hpp binomial_dist.hpp negative_binomial_dist.hpp \
		     hypergeometric_dist.hpp geometric_dist.hpp \
		     poisson_dist.hpp zero_truncated_poisson_dist.hpp discrete_dist.hpp fast_discrete_dist.hpp

libtrng4_la_LDFLAGS = -version-info 19:0:0
