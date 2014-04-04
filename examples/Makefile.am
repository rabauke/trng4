AUTOMAKE_OPTIONS = gnu

AM_CPPFLAGS = -I$(top_builddir)
AM_LDFLAGS = -L../src/.libs -ltrng4 -static
AM_MPILDFLAGS = -L../src/.libs -ltrng4
AM_CXXFLAGS = @BOOST_CXXFLAGS@ @CXXFLAGS@

EXTRA_PROGRAMS = hello_world sample_output time pi pi_block_mpi pi_block_openmp pi_block_tbb pi_leap_mpi pi_leap_openmp stl_container cpp11 distributions plausibility_test Ising_model correlated_normal_dist discrete_dist discrete_dist_c_style bernoulli_dist

EXTRA_DIST = pi_block_cuda.cu pi_leap_cuda.cu

MPICXXCOMPILE = $(MPICXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
MPICXXLD = $(MPICXX)
MPICXXLINK = $(MPICXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) -o $@ 

hello_world_SOURCES = hello_world.cc

sample_output_SOURCES = sample_output.cc

time_SOURCES = time.cc

pi_SOURCES = pi.cc

pi_block_mpi_SOURCES = pi_block_mpi.cc

pi_block_openmp_SOURCES = pi_block_openmp.cc

pi_block_tbb_SOURCES = pi_block_tbb.cc

pi_leap_mpi_SOURCES = pi_leap_mpi.cc

pi_leap_openmp_SOURCES = pi_leap_openmp.cc

stl_container_SOURCES = stl_container.cc

cpp11_SOURCES = cpp11.cc

distributions_SOURCES = distributions.cc

plausibility_test_SOURCES = plausibility_test.cc

Ising_model_SOURCES = Ising_model.cc

correlated_normal_dist_SOURCES = correlated_normal_dist.cc

discrete_dist_SOURCES = discrete_dist.cc

discrete_dist_c_style_SOURCES = discrete_dist_c_style.cc

bernoulli_dist_SOURCES = bernoulli_dist.cc

$(pi_block_mpi_OBJECTS): $(pi_block_mpi_SOURCES)
	$(MPICXXCOMPILE) -c pi_block_mpi.cc

pi_block_mpi$(EXEEXT): $(pi_block_mpi_OBJECTS) $(pi_block_mpi_DEPENDENCIES)
	$(MPICXXLINK) $(pi_block_mpi_OBJECTS) $(AM_MPILDFLAGS) 

$(pi_block_openmp_OBJECTS): $(pi_block_openmp_SOURCES)
	$(CXXCOMPILE) $(OPENMP_CXXFLAGS) -c pi_block_openmp.cc

pi_block_openmp$(EXEEXT): $(pi_block_openmp_OBJECTS) $(pi_block_openmp_DEPENDENCIES)
	$(CXXLINK) $(OPENMP_CXXFLAGS) $(pi_block_openmp_LDFLAGS) $(pi_block_openmp_OBJECTS) $(pi_block_openmp_LDADD) $(LIBS)

$(pi_block_tbb_OBJECTS): $(pi_block_tbb_SOURCES)
	$(CXXCOMPILE) -c pi_block_tbb.cc

pi_block_tbb$(EXEEXT): $(pi_block_tbb_OBJECTS) $(pi_block_tbb_DEPENDENCIES)
	$(CXXLINK) $(pi_block_tbb_LDFLAGS) $(pi_block_tbb_OBJECTS) $(pi_block_tbb_LDADD) $(LIB_TBB) $(LIBS)

$(pi_leap_mpi_OBJECTS): $(pi_leap_mpi_SOURCES)
	$(MPICXXCOMPILE) -c pi_leap_mpi.cc

pi_leap_mpi$(EXEEXT): $(pi_leap_mpi_OBJECTS) $(pi_leap_mpi_DEPENDENCIES)
	$(MPICXXLINK) $(pi_leap_mpi_OBJECTS) $(AM_MPILDFLAGS)

$(pi_leap_openmp_OBJECTS): $(pi_leap_openmp_SOURCES)
	$(CXXCOMPILE) $(OPENMP_CXXFLAGS) -c pi_leap_openmp.cc

pi_leap_openmp$(EXEEXT): $(pi_leap_openmp_OBJECTS) $(pi_leap_openmp_DEPENDENCIES)
	$(CXXLINK) $(OPENMP_CXXFLAGS) $(pi_leap_openmp_LDFLAGS) $(pi_leap_openmp_OBJECTS) $(pi_leap_openmp_LDADD) $(LIBS)


examples: $(EXTRA_PROGRAMS)

clean-local:
	-rm -f $(EXTRA_PROGRAMS)
