# -------------------------------
# Compilers
# -------------------------------
# C++ compiler
CXX = g++

# Fortran compiler
FC  = gfortran


# Detect operating system (Linux / Darwin / Windows)
UNAME_S := $(shell uname -s)


# -------------------------------
# OS-specific settings
# -------------------------------
# These control:
# - library extension (.so / .dll / .dylib)
# - position-independent code (needed for shared libs)
# - linker flags
# - remove command

ifeq ($(UNAME_S),Linux)
    EXT = so            # shared library extension
    PIC = -fPIC         # required for shared libraries
    SHARED = -shared    # build shared library
    RM = rm -f
endif

ifeq ($(UNAME_S),Darwin)
    EXT = dylib         # macOS shared library extension
    PIC = -fPIC
    SHARED = -shared
    RM = rm -f
endif

# Windows (Rtools / MinGW)
ifeq ($(OS),Windows_NT)
    EXT = dll
    PIC =               # not needed on Windows
    SHARED = -shared
    RM = rm -f          # works under Rtools
endif


# -------------------------------
# Compiler flags
# -------------------------------
# -O3        : maximum optimization
# -fopenmp   : enable OpenMP parallelism
# $(PIC)     : position-independent code (Linux/macOS)
# -Iinclude  : look for header files in ./include

CXXFLAGS = -O3 -fopenmp $(PIC) -Iinclude
FFLAGS   = -O3 -fopenmp $(PIC)


# -------------------------------
# Targets (what we build)
# -------------------------------
# Combined library (C++ + Fortran)
LIB     = libqgtools.$(EXT)

# C++-only library
CPP_LIB = libqgc.$(EXT)

# Fortran-only library
F_LIB   = libqgf.$(EXT)

# Object files (compiled but not yet linked)
OBJS = grsbed_c.o grsbed_f.o


# Default target: build everything
all: $(CPP_LIB) $(F_LIB) $(LIB)


# -------------------------------
# Combined library
# -------------------------------
# Build one shared library containing BOTH:
#   - C++ code
#   - Fortran code
#
# IMPORTANT:
# We link with gfortran (FC), not g++, because:
#   - gfortran automatically links Fortran runtime
#   - we manually add C++ runtime via -lstdc++
#
# This ensures both:
#   - libgfortran  (Fortran runtime)
#   - libstdc++    (C++ runtime)
# are available

$(LIB): $(OBJS)
	$(FC) $(SHARED) $(OBJS) -o $(LIB) -fopenmp -lstdc++


# -------------------------------
# Object file compilation
# -------------------------------
# Compile C++ source → object file
# -c = compile only (no linking yet)

grsbed_c.o: grsbed_c.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile Fortran source → object file
grsbed_f.o: grsbed_f.f90
	$(FC) $(FFLAGS) -c $< -o $@


# -------------------------------
# Separate libraries (for demo purposes)
# -------------------------------
# These illustrate that we can build each backend independently

# Build C++-only shared library
$(CPP_LIB): grsbed_c.cpp
	$(CXX) $(CXXFLAGS) $(SHARED) $< -o $@

# Build Fortran-only shared library
$(F_LIB): grsbed_f.f90
	$(FC) $(FFLAGS) $(SHARED) $< -o $@


# -------------------------------
# Clean build artifacts
# -------------------------------
# Remove all compiled files

clean:
	$(RM) *.o *.so *.dll *.dylib
