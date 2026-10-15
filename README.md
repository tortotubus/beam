# Basilisk Beam Library

[![Unit Tests](https://github.com/tortotubus/beam/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/tortotubus/beam/actions/workflows/unit-tests.yml)


## Dependencies
For RPM-based distros (e.g. Redhat, Fedora):
```bash
sudo dnf install -y cmake gcc gcc-c++ git gcc-gfortran make gawk gnuplot doxygen gtest-devel openmpi-devel eigen3-devel hdf5-devel hdf5-openmpi-devel
```

For apt-based distros (e.g. Ubuntu):
```bash
sudo apt-get install -y  cmake gcc g++ git make gawk gnuplot doxygen libgtest-dev libopenmpi-dev libeigen3-dev libhdf5-dev  
```

## Building
From the repository root:
```bash
cmake -S . --preset Default 
cmake --build --preset Default --parallel  
```
 
