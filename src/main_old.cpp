#include <beam/Euler/StaticEulerBeam.hpp>
#include <beam/FiniteElement/Mesh/MeshAdapters/VTKAdapter.hpp>
#include <beam/FiniteElement/Mesh/MeshAdapters/VTKWindow.hpp>
#include <beam/FiniteElement/Mesh/MeshReaders/PLYReader.hpp>

#include "mfem.hpp"

using namespace beam;

int test_ply_reader() {
  // Load the Stanford Bunny (ascii PLY) into a triangle mesh:
  auto bunny = beam::PLYReader::load(
      "/home/colive/Downloads/bunny/reconstruction/bun_zipper.ply");

  std::cout << "Loaded bunny: " << bunny.size_nodes() << " vertices, "
            << bunny.size_elements() << " triangles.\n";

  // Create and show visualization window
  VTKWindow window(beam::toVTKGrid(bunny));
  window.set_background(1.0, 1.0, 1.0);  // white
  window.show();

  return 0;
}

int test_static_beam() {
  const double load = 1.;
  const double length = 1.; // Length of beam
  const double EI = 1.;     // Bending stiffness
  const size_t N = 20;      // Number of elements

  beam::StaticEulerBeam beam(load, EI, length, N);

  beam.get_mesh().plot_gnuplot();
  VTKWindow window(beam::toVTKGrid(beam.get_mesh()));
  window.set_background(1.0, 1.0, 1.0);  // white
  window.show();

  return 0;
}

int main() { return test_static_beam(); }