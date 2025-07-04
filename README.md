# Basilisk Beam Library

## Overview
This C++ library provides tools for beam mechanics calculations, including:
- Matrix and Vector operations
- Euler beam solvers
- Finite element implementations

## Dependencies
For RPM-based distros:
```bash
sudo dnf install -y cmake gcc gcc-c++ gcc-gfortran make gawk gnuplot doxygen gtest-devel openmpi-devel eigen3-devel hdf5-devel hdf5-openmpi-devel
```
On sockeye.arc.ubc.ca:
```bash
module load gcc openmpi hdf5 cmake eigen
```

## Building
```bash
mkdir build && cd build
cmake ..
make
```

## Documentation
Generate documentation using:
```bash
make doc
```
Then open `docs/html/index.html` in your browser.
