# -------------------------------
# Compilers
# -------------------------------
CXX = g++
FC  = gfortran
CC  = gcc

UNAME_S := $(shell uname -s)

# -------------------------------
# OS-specific settings
# -------------------------------
ifeq ($(UNAME_S),Linux)
    EXT = so
    PIC = -fPIC
    SHARED = -shared
    RM = rm -f
endif

ifeq ($(UNAME_S),Darwin)
    EXT = dylib
    PIC = -fPIC
    SHARED = -dynamiclib
    RM = rm -f
endif

ifeq ($(OS),Windows_NT)
    EXT = dll
    PIC =
    SHARED = -shared
    RM = rm -f
endif

# -------------------------------
# R config
# -------------------------------
R_CPPFLAGS := $(shell R CMD config --cppflags)
R_LDFLAGS  := $(shell R CMD config --ldflags)

# -------------------------------
# Flags
# -------------------------------
CXXFLAGS = -O3 -fopenmp $(PIC) $(R_CPPFLAGS)
FFLAGS   = -O3 -fopenmp $(PIC)
CFLAGS   = $(R_CPPFLAGS)

# -------------------------------
# Targets
# -------------------------------
LIB = libqgtools.$(EXT)

# IMPORTANT: include wrapper + init
OBJS = \
  mtgrsbed_core.o \
  mtgrsbed_core_c.o \
  mtgrsbed_core_f.o \
  mtgrsbed_call.o \
  libqgtools_init.o


all: $(LIB)

# -------------------------------
# Link shared library (WITH R!)
# -------------------------------
$(LIB): $(OBJS)
	$(FC) $(SHARED) $(OBJS) -o $(LIB) \
	    -fopenmp -lstdc++


# -------------------------------
# Compilation rules
# -------------------------------
mtgrsbed_core.o: mtgrsbed_core.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

mtgrsbed_core_c.o: mtgrsbed_core_c.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

mtgrsbed_core_f.o: mtgrsbed_core_f.f90
	$(FC) $(FFLAGS) -c $< -o $@

mtgrsbed_call.o: mtgrsbed_call.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

libqgtools_init.o: libqgtools_init.c
	$(CC) $(CFLAGS) -c $< -o $@

# -------------------------------
# Clean
# -------------------------------
clean:
	$(RM) *.o *.so *.dll *.dylib
