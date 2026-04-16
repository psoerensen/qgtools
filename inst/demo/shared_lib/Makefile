CXX = g++
FC  = gfortran

CXXFLAGS = -O3 -fopenmp -fPIC
FFLAGS   = -O3 -fopenmp -fPIC

all: libgrsbed_cpp.so libgrsbed_fortran.so

libgrsbed_cpp.so: grsbed.cpp
	$(CXX) $(CXXFLAGS) -shared grsbed.cpp -o libgrsbed_cpp.so

libgrsbed_fortran.so: grsbed_fortran.f90
	$(FC) $(FFLAGS) -shared grsbed_fortran.f90 -o libgrsbed_fortran.so

clean:
	rm -f *.o *.so
