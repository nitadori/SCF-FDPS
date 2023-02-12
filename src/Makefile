PS_PATH = ../FDPS/src/
PS_INC = -I$(PS_PATH)

CXX = g++
MPICXX = mpicxx
CXXFLAGS = -std=c++17 -O3 -march=native
CXXFLAGS += -Wall
#CLIBS = -ltcmalloc
#CXXFLAGS += -ffast-math
#CXXFLAGS += -funroll-loops
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
PS_DEFS = -DPARTICLE_SIMULATOR_MPI_PARALLEL
PS_DEFS += -DQUADRUPOLE

scfdhmpi_v0.out:scfdhmpi_v0.o scfio.o scfmpi.o scfkernel.o scfmath.o scfstep.o scfcoef.s fp_check.s gravity_kernel.s gravity_spline.s
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(CLIBS)

scfdhmpi_v0.o:scfdhmpi_v0.cpp user_defined.hpp
	time $(MPICXX) $(PS_INC) $(CXXFLAGS) $(PS_DEFS) -c $<

scfio.o:scfio.cpp
	$(MPICXX) $(CXXFLAGS) -c $^

scfmpi.o:scfmpi.cpp
	$(MPICXX) $(CXXFLAGS) -c $^

scfkernel.o:scfkernel.cpp
	$(CXX) $(CXXFLAGS)  -c $^

scfmath.o:scfmath.cpp
	$(CXX) $(CXXFLAGS) -c $^

scfstep.o:scfstep.cpp
	$(CXX) $(CXXFLAGS) -c $^

difftest.out: difftest.s scfcoef.s scfkernel.o scfmath.o scfstep.o scfio.o scfmpi.o
	$(MPICXX) $(CXXFLAGS) -Ofast -o $@ $^

difftest.s: difftest.cpp
	$(MPICXX) $(CXXFLAGS) -Ofast -S $^

scfcoef.s:scfcoef.cpp scfcoef-inc.h
	$(CXX) $(CXXFLAGS) -Ofast -S $<

fp_check.s: fp_check.c
	gcc -O1 -march=native -S $<

gravity_kernel.s: gravity_kernel.cpp
	$(CXX) $(CXXFLAGS) $(PS_INC) -S $^

gravity_spline.s: gravity_spline.cpp
	$(CXX) $(CXXFLAGS) $(PS_INC) -S $^

diffaccp.out: diffaccp.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

benchmark.out: benchmark.cpp gravity_kernel.s gravity_spline.s
	$(CXX) $(CXXFLAGS) $(PS_INC) -o $@ $^

b2ah.out: b2ah.c
	gcc -O -Wall -o $@ $<

b2ad.out: b2ad.c
	gcc -O -Wall -o $@ $<

all: scfdhmpi_v0.out difftest.out diffaccp.out benchmark.out b2ah.out b2ad.out

clean:
	rm -f *.o *.s *.out

version:
	which $(CXX)
	$(CXX) --version
	which $(MPICXX)
	$(MPICXX) --version
