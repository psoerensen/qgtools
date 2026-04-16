This folder demonstrates how to build standalone shared libraries
from qgtools computational kernels.

Contents:
- C++ OpenMP implementation
- Fortran OpenMP implementation
- R example calling both

Build:
  make

Run:
  Rscript test.R

Purpose:
- illustrate shared library design
- show C++ ↔ Fortran ↔ R interoperability
- demonstrate OpenMP parallel execution
