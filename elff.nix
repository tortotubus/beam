{
  lib,
  stdenv,
  cmake,
  pkg-config,
  gnumake,
  gcc,
  openmpi,
  gtest,
  eigen,
  hdf5-mpi,
  pacific,
  criterion,
  makeWrapper,
}:

stdenv.mkDerivation {
  pname = "elff";
  version = "0.0.1";

  # Build the current source tree (this repo checkout)
  src = ./.;

  nativeBuildInputs = [
    gcc
    cmake
    gnumake
    pkg-config
    criterion
    pacific
    makeWrapper
  ];

  buildInputs = [
    gtest
    eigen
    (hdf5-mpi.override { mpi = openmpi; })
    openmpi
  ];

  propagatedBuildInputs = [
    gtest
    eigen
    (hdf5-mpi.override { mpi = openmpi; })
    openmpi
  ];

  cmakeFlags = [
    "-DCMAKE_BUILD_TYPE=Release"
    "-DBUILD_BASILISK_EXAMPLES=OFF"
  ];

  meta = with lib; {
    description = "ELFF: ELastic Fibers in Fluid";
    homepage = "https://github.com/tortotubus/PacIFiC";
    platforms = platforms.linux;
    license = licenses.mit;
  };
}
