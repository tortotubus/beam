#define MPI_AUTO_BOUNDARY 0

// #include "grid/quadtree.h"
#include "grid/multigrid.h"

#include "library/ibm/IBKernels.h"
#include "library/ibm/IBMeshManager.h"

#if TREE
#include "library/ibm/IBAdapt.h"
#endif

#include "library/io/output-vtk.h"

const int maxlevel = 6;
const int minlevel = 3;
int ibmlevel = 10;

scalar suppf[];
scalar weightf[];
scalar levelf[];
scalar pidf[];
scalar gradientf[];

scalar ibmf[];

IBscalar ibgradientf;
IBscalar iblevelf;
IBscalar ibnodeidf;

coord
circle(int n, int N, coord centre, double radius)
{
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}

#if _MPI
static void
print_exchange_summary(void)
{

  if (pid() == 0){
    printf("GHOSTS: %d \n", GHOSTS);
    printf("BGHOSTS: %d \n", BGHOSTS);
  }

  for (int rank = 0; rank < npe(); rank++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == pid()) {
      printf("[recv %d] {", pid());
      bool first = true;
      for (int peer = 0; peer < npe(); peer++) {
        if (peer == pid())
          continue;
        int size = ibmm.rcv_boundary[peer].nodes.size;
        if (size > 0) {
          for (int i = 0; i < size; i++) {
            if (first) {
              printf("%d", peer);
              first = false;
            } else {
              printf(" %d", peer);
            }
          }
        }
      }
      printf("}\n");

      printf("[snd %d] {", pid());
      first = true;
      for (int peer = 0; peer < npe(); peer++) {
        if (peer == pid())
          continue;
        int size = ibmm.snd_boundary[peer].nodes.size;
        if (size > 0) {
          for (int i = 0; i < size; i++) {
            if (first) {
              printf("%d", peer);
              first = false;
            } else {
              printf(" %d", peer);
            }
          }
        }
      }
      printf("}\n");
      fflush(stdout);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

int
main()
{
  new_ibscalar(ibgradientf);
  new_ibscalar(iblevelf);
  new_ibscalar(ibnodeidf);

  X0 = -4, Y0 = -4, L0 = 8;
  // periodic(right);
  // periodic(top);
#if TREE
  init_grid(1 << minlevel);
#else
  init_grid(1 << ibmlevel);
  ibmlevel = depth();
#endif

  // printf("%d %d \n", ibmlevel, depth());

  // Create mesh manager object with 0 meshes
  ibmeshmanager_init(0);

  // Add meshes
  int new_id = ibmeshmanager_add_mesh();

  // Add a single node to the first mesh
  const int N_circ = 100;
  const double r_circ = .15;
  const coord cen_circ = { 0.0 };

  ibmeshmanager_add_nodes(new_id, N_circ);

  foreach_ibnode()
  {
    node->depth = depth();
  }

  foreach_ibnode_per_ibmesh()
  {
    node->pos = circle(node_id, N_circ, cen_circ, r_circ);
  }

  // Refine and coarsen around nodes
#if TREE
  adapt_wavelet_ibm(NULL, NULL, maxlevel, minlevel);
#endif

#if _MPI
  ibmeshmanager_update_pid();
#endif

  foreach () {
    levelf[] = point.level;
    suppf[] = 0.;
    weightf[] = 0.;
    pidf[] = pid();
#if dimension == 1
    gradientf[] = x;
#elif dimension == 2
    gradientf[] = x + y;
#else
    gradientf[] = x + y + z;
#endif
  }

  foreach_ibnode()
  {
    ibval(ibnodeidf) = node_id;
  }

  foreach_ibnode_per_ibmesh()
  {
    peskin_cosine_kernel_spread_dimensionless(node)
    {
      suppf[] = pid() + 1;
      weightf[] += weight;
    }
  }

#if _MPI
  gradientf.dirty = true;
  levelf.dirty = true;
  boundary((scalar*){ gradientf, levelf });
#endif

  foreach_ibnode()
  {
    peskin_cosine_kernel_gather_dimensionless(node)
    {
      ibval(ibgradientf) += weight * gradientf[];
      ibval(iblevelf) += weight * levelf[];
    }
  }

#if _MPI
  print_exchange_summary();
  ibmeshmanager_boundary();
#endif

  foreach_ibnode()
  {
    int cell_count = 0;
    peskin_cosine_kernel_spread_dimensionless(node)
    {
      ibmf[] += ibval(ibgradientf) * weight;
      cell_count++;
    }
  }
 

  // Output
  {
    int i = 0;
    double t = 0;
#if TREE
    output_hdf_htg();
#else
    output_hdf_imagedata();
#endif
    output_hdf_pd();
  }

  ibmeshmanager_free();

  return 0;
}
