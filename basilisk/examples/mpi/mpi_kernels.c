#include "grid/quadtree.h"
// #include "grid/multigrid.h"

#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBKernels.h"

#if TREE
  #include "library/ibm/IBAdapt.h"
#endif 

#include "library/io/output-vtk.h"

const int maxlevel = 6;
const int minlevel = 3;
const int ibmlevel = 10;

scalar suppf[];
scalar weightf[];
scalar levelf[];
scalar pidf[];
scalar gradientf[];

coord
circle(int n, int N, coord centre, double radius)
{
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}


static void mpi_write_neighbors (const char * basename)
{
  if (npe() <= 1)
    return;

  // Ensure the lists are up-to-date (recomputed in mpi_boundary_update_buffers()) :contentReference[oaicite:1]{index=1}
  mpi_boundary_update_buffers();

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  if (!m || !m->send || !m->receive)
    return;

  char fname[256];
  snprintf(fname, sizeof(fname), "%s-%d.txt", basename, pid());
  FILE * fp = fopen(fname, "w");
  if (!fp) return;

  int ns = m->send->len/sizeof(int);
  int nr = m->receive->len/sizeof(int);
  int * snd = (int *) m->send->p;
  int * rcv = (int *) m->receive->p;

  fprintf(fp, "pid %d / %d\n", pid(), npe());
  fprintf(fp, "send neighbors (%d):", ns);
  for (int i = 0; i < ns; i++) fprintf(fp, " %d", snd[i]);
  fprintf(fp, "\n");

  fprintf(fp, "recv neighbors (%d):", nr);
  for (int i = 0; i < nr; i++) fprintf(fp, " %d", rcv[i]);
  fprintf(fp, "\n");

  fclose(fp);
}

int
main()
{
#if TREE
  init_grid(1 << minlevel);
#else 
  init_grid(1 << ibmlevel);
#endif

  X0 = -4, Y0 = -4, L0 = 8;

  periodic(right);
  periodic(top);

  // Create mesh manager object with 0 meshes
  ibmeshmanager_init(0);

  // Add meshes
  int new_id = ibmeshmanager_add_mesh();

  // Add a single node to the first mesh
  const int N_circ = 100;
  const double r_circ = .15;
  const coord cen_circ = { 0.0 };

  ibmeshmanager_add_nodes(new_id, N_circ);

  foreach_ibmesh()
  {
    mesh->refinement_level = ibmlevel;
  }

  foreach_ibnode_per_ibmesh()
  {
    node->lagpos = circle(node_id, N_circ, cen_circ, r_circ);
  }


  // Refine and coarsen around nodes 
#if TREE
  adapt_wavelet_ibm(NULL, NULL, maxlevel, minlevel);
#endif

  ibmeshmanager_init_stencil_caches();

  foreach() {
    levelf[] = point.level;
    suppf[] = 0.;
    weightf[] = 0.;
    pidf[] = pid();
  #if dimension == 1
    gradientf[] = x;
  #elif dimension == 2 
    gradientf[] = x+y;
  #else
    gradientf[] = x+y+z;
  #endif
  }

  foreach_ibnode() {
    ibnode_update_pid(node);
  }

  foreach_ibnode() {
    peskin_cosine_kernel_spread_dimensionless(node) {
      suppf[] = pid() + 1;
      weightf[] += weight;
    }
  }

#if _MPI
  gradientf.dirty = true;
  levelf.dirty = true;
  boundary((scalar*){gradientf, levelf});
#endif

  foreach_ibnode() {
    peskin_cosine_kernel_gather_dimensionless(node) {
      node->force.x += weight * gradientf[];
      node->force.y += weight * levelf[];
    }
  }

  // Output
  {
    int i = 0;
    double t = 0;
    output_hdf_htg();
    output_hdf_pd();
  }

  if (pid() == 5) {
    foreach_ibnode() {
      printf("node pid: %d\n", node->pid);
    }
  }

  ibmeshmanager_free();

  return 0;
}
