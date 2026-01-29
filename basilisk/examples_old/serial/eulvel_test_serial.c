#include "library/ibm/IBcentered.h"
#include "library/ibm/IBoutput.h"
#include "library/output_htg.h"

int main() {

  L0 = 8;
  init_grid(1 << 9);
  origin(0,0);

  foreach() {
    p[] = pid();
    u.x[] = 1;
    u.y[] = 1;
  }

  int n_meshes = 1;
  ib_mesh_manager_init(n_meshes);
 
  double b_length = 1;  
  int b_nodes = 65;
  double b_ds = b_length / b_nodes; 
  double b_theta = 0.1 * pi;
 
  ib_mesh_init_nn(ibmesh(0),b_nodes);

  foreach_ibnode(ibmesh(0))
  {
    double radius = (node_index * b_ds);
    node->lagpos.x = L0/2 -0.5 + radius * cos(b_theta);
    node->lagpos.y = L0/2 -0.5 + radius * sin(b_theta);
    node->weight = 1;
  }

  ibmesh(0)->isactive = true;
  ib_mesh_manager_generate_stencil_cache();
  ib_mesh_manager_interpolate_eulerian_velocities();

  output_hdf_htg({p}, {u}, "interpolation_test_serial", 0, 0);
  output_ibnodes("interpolation_test_serial",0,0);
}