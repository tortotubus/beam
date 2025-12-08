#define DEBUG 0

int maxlevel = 9;
int Nf = 1 << 9;

#include "grid/quadtree.h"

#include "interface/ibm/beam.h"
#include "library/ibm/uzawa2.h"

int
main()
{
  L0 = 8;
  origin(-2, -4.);
  init_grid(Nf);
  DT = 0.0005;

  int n_meshes = 2;
  ib_mesh_manager_init(n_meshes);

  double f_dx = L0/Nf;
  // double b_dx = f_dx * 1.683;

  double b_length = 1;
  double b_EI = 0.0001;
  double b_mu = 1.5;
  int b_nodes = 15;
  double b_ds = b_length/b_nodes;
  double b_r = 1e4;
  double b_theta = 0.1*pi;

  ib_beam_t beam1 = ib_beam_new_theta(2, b_EI, b_mu, b_nodes, b_r, b_theta);
  ib_mesh_model(ibmesh(0), beam1);
  ibmesh(0)->pid = 0;

  ib_beam_t beam2 = ib_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta);
  ib_mesh_model(ibmesh(1), beam2);
  ibmesh(1)->pid = 1;

  foreach_ibmesh() {
    foreach_ibnode(mesh) {
      node->gravity.x = 1;
    }
  }

  // ib_mesh_manager_set_owners();

  for (int i = 0; i < 50; i++) {
    ib_mesh_manager_advance_lagrangian_mesh();
  }

  foreach_ibmesh() {
    vertex_t centroid = ib_mesh_get_centroid(mesh);
    printf(
      "[pid %d : mesh %d]: %f,%f,%f\n", pid(), mesh_index, centroid.x, centroid.y, centroid.z);
  }
}