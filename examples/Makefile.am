AUTOMAKE_OPTIONS = gnu

INCLUDES = -I$(top_builddir)

LDADD = ../src/.libs/libtrng4.a

EXTRA_PROGRAMS = hello_world test time pi pi_block_mpi pi_block_openmp pi_leap_mpi pi_leap_openmp stl_container distributions implementation_test Ising_model correlated_normal_dist
EXTRA_OBJECTS = $(addsuffix .o, $(EXTRA_PROGRAMS))

MPICXXCOMPILE = $(MPICXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
MPICXXLD = $(MPICXX)
MPICXXLINK = $(MPICXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@

hello_world_SOURCES = hello_world.cc

test_SOURCES = test.cc

time_SOURCES = time.cc

pi_SOURCES = pi.cc

pi_block_mpi_SOURCES = pi_block_mpi.cc

pi_block_openmp_SOURCES = pi_block_openmp.cc

pi_leap_mpi_SOURCES = pi_leap_mpi.cc

pi_leap_openmp_SOURCES = pi_leap_openmp.cc

stl_container_SOURCES = stl_container.cc

distributions_SOURCES = distributions.cc

implementation_test_SOURCES = implementation_test.cc

Ising_model_SOURCES = Ising_model.cc

correlated_normal_dist_SOURCES = correlated_normal_dist.cc


$(pi_block_mpi_OBJECTS): $(pi_block_mpi_SOURCES)
	$(MPICXXCOMPILE) -c pi_block_mpi.cc

pi_block_mpi$(EXEEXT): $(pi_block_mpi_OBJECTS) $(pi_block_mpi_DEPENDENCIES)
	$(MPICXXLINK) $(pi_block_mpi_LDFLAGS) $(pi_block_mpi_OBJECTS) $(pi_block_mpi_LDADD) $(LIBS)


$(pi_block_openmp_OBJECTS): $(pi_block_openmp_SOURCES)
	$(CXXCOMPILE) $(OPENMP_CXXFLAGS) -c pi_block_openmp.cc

pi_block_openmp$(EXEEXT): $(pi_block_openmp_OBJECTS) $(pi_block_openmp_DEPENDENCIES)
	$(CXXLINK) $(OPENMP_CXXFLAGS) $(pi_block_openmp_LDFLAGS) $(pi_block_openmp_OBJECTS) $(pi_block_openmp_LDADD) $(LIBS)


$(pi_leap_mpi_OBJECTS): $(pi_leap_mpi_SOURCES)
	$(MPICXXCOMPILE) -c pi_leap_mpi.cc

pi_leap_mpi$(EXEEXT): $(pi_leap_mpi_OBJECTS) $(pi_leap_mpi_DEPENDENCIES)
	$(MPICXXLINK) $(pi_leap_mpi_LDFLAGS) $(pi_leap_mpi_OBJECTS) $(pi_leap_mpi_LDADD) $(LIBS)

$(pi_leap_openmp_OBJECTS): $(pi_leap_openmp_SOURCES)
	$(CXXCOMPILE) $(OPENMP_CXXFLAGS) -c pi_leap_openmp.cc

pi_leap_openmp$(EXEEXT): $(pi_leap_openmp_OBJECTS) $(pi_leap_openmp_DEPENDENCIES)
	$(CXXLINK) $(OPENMP_CXXFLAGS) $(pi_leap_openmp_LDFLAGS) $(pi_leap_openmp_OBJECTS) $(pi_leap_openmp_LDADD) $(LIBS)


examples: $(EXTRA_PROGRAMS)

clean:
	-rm -f $(EXTRA_PROGRAMS) $(EXTRA_OBJECTS)
