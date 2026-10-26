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
  libxml2,
  ncurses,
  basiliskSrc,
  criterion,
  gcovr,
  makeWrapper,
}:

stdenv.mkDerivation {
  pname = "elff";
  version = "0.0.1";

  # Build the current source tree (this repo checkout)
  src = ./.;

  nativeBuildInputs = [
    gcc
    gcovr
    cmake
    gnumake
    pkg-config
    criterion
    makeWrapper
  ];

  buildInputs = [
    gtest
    eigen
    (hdf5-mpi.override { mpi = openmpi; })
    libxml2
    ncurses
    openmpi
  ];

  propagatedBuildInputs = [
    gtest
    eigen
    (hdf5-mpi.override { mpi = openmpi; })
    libxml2
    openmpi
  ];

  # Flake Git snapshots do not include submodule worktrees. Populate the
  # submodule path from the separately pinned source before configuring CMake.
  postPatch = ''
    mkdir -p third_party/basilisk/basilisk_submodule
    cp -R ${basiliskSrc}/. third_party/basilisk/basilisk_submodule/
  '';

  cmakeFlags = [
    "-DCMAKE_BUILD_TYPE=Release"
    "-DELFF_USE_SUBMODULES=ON"
    "-DBUILD_BASILISK_EXAMPLES=ON"
  ];

  shellHook = ''
    export GCOV="gcov -i"
  '';

  meta = with lib; {
    description = "ELFF: ELastic Fibers in Fluid";
    homepage = "https://github.com/tortotubus/PacIFiC";
    platforms = platforms.linux;
    license = licenses.mit;
  };
}
