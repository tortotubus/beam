@page build Building

## Dependencies

This library has a handfull of dependencies which are common enough to be availible either through your repository package manager (if building on a desktop computer) or otherwise through the module system in your supercomputing environment. The hard-dependencies for the library itself are
  - ``gcc`` (or any other modern C++ compiler)
  - ``openmpi``
  - ``cmake``
  - ``eigen3``
  - ``hdf5``

Ideally, optional dependencies may include 

  - Any linear algebra backend, for ``eigen3``
  - Google Test, to run unit tests
  - Doxygen, for building this documentation
  - gnuplot, for plotting some results from the unit tests
  
In addition, a lightly modified and forked copy of [basilisk](http://basilisk.fr/) can be used in this repository and is provided with an immersed boundary coupling to some elastic models in this library. The main modifications are replacement of the ``make``-based build system with a ``cmake`` build system, as well as including ``hdf5`` libraries so that the octree fluid grid may be saved to disk as a [HyperTreeGrid](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#hypertreegrid) which is supported in Paraview >6.0. If you choose to build it, you will need to install its own depencies, which include 

  - ``gcc``
  - ``openmpi``
  - ``cmake``
  - ``hdf5``

The main target that basilisk provides is ``qcc``, a meta-compiler which compiles "Basilisk C" to C code, which is then compiled using any C99 compiler. Normally, ``qcc`` makes this transparent to the user, unless the flag ``-source`` is used, in which case ``qcc`` outputs a portable C99 source file rather than a binary. The CMake build system provided in this library will fetch the basilisk source from [here](https://github.com/tortotubus/beam-basilisk), build ``qcc``, and use this to compile the example basilisk programs in the ``basilisk/examples`` directory, automatically linking to the beam library which exposes a small C interface defined in ``basilisk/interface``.

## Building


## Installing

