#include "library/df/kernels.h"
#include "navier-stokes/centered.h"

#if dimension == 2
#define STENCIL_SIZE 25
#else // dimension == 3
#define STENCIL_SIZE 125
#endif

/**
 * @brief Contains information for a single node of the immersed boundary method
 */
typedef struct
{
  coord lagpos;  /** Lagrangian position */
  coord force;     /** Force density on the node */
  coord lagvel;    /** Lagrangian velocity */
  coord eulvel;    /** Interpolated Eulerian velocity */
  //double residual; /** Residual between Lagrangian and Eulerian velocity */
  Cache stencil; /** Stencil is used to cache all grid cells in the neighborhood
                  of a given node, so that they do not need to be located
                  multilpe times */
} IBNode;

/**
 * @brief Initialize an IBNode struct
 */
void
ib_node_init(IBNode* n)
{
  n->lagpos.x = 0.;
  n->lagpos.y = 0.;

  n->force.x = 0.;
  n->force.y = 0.;

  n->lagvel.x = 0.;
  n->lagvel.y = 0.;

  n->eulvel.x = 0.;
  n->eulvel.y = 0.;

  n->stencil.n = 0.;
  n->stencil.nm = STENCIL_SIZE;
  n->stencil.p = (Index*)malloc(STENCIL_SIZE * sizeof(Index));
}

/**
 * @brief Safely free the members of an IBNode struct
 */
void
ib_node_free(IBNode* n)
{
  free(n->stencil.p);
}

/**
 * @brief Represent an immersed boundary mesh as a collection of IBNodes
 */
typedef struct
{
  IBNode* nodes; /** Heap-allocated nodes */
  int nn;        /** The number of heap-allocated nodes */
  bool isactive; /** Tells the program if the mesh is active */
} IBMesh;

/**
 * @brief Safely free an IBMesh
 *
 * @param mesh Pointer to the mesh
 */
void
ib_mesh_free(IBMesh* mesh)
{
  if (mesh->nodes) {
    for (int ni = 0; ni < mesh->nn; ni++) {
      ib_node_free(&mesh->nodes[ni]);
    }
    free(mesh->nodes);
  }
}

/**
 * @brief Initialize an IBMesh
 *
 * @param mesh Pointer to the mesh
 * @param nn Number of nodes in the mesh
 */
void
ib_mesh_init(IBMesh* mesh, int nn)
{
  mesh->nodes = (IBNode*)calloc(nn, sizeof(IBNode));
  mesh->nn = nn;
  mesh->isactive = true;
  for (int ni = 0; ni < nn; ni++) {
    ib_node_init(&mesh->nodes[ni]);
  }
}

/**
 * @brief Management struct to hold all IBMeshes in the simulation
 */

typedef struct
{
  IBMesh* meshes;
  int nm;
} IBMeshManager;

#define IBMESH(i) (ib_mesh_manager.meshes[i])

IBMeshManager ib_mesh_manager;

/**
 * @brief Initialize the immersed boundary mesh manager
 *
 * @param nm The number of meshes you plan to have
 */
void
ib_mesh_manager_init(int nm)
{
  if (ib_mesh_manager.meshes == NULL) {
    ib_mesh_manager.nm = nm;
    ib_mesh_manager.meshes = (IBMesh*)calloc(nm, sizeof(IBMesh));
  } else {
    return;
  }
}

/**
 * @brief Free all members in the immersed boundary mesh manager, including the
 * array of meshes, their nodes, and stencil nodes.
 */
void
ib_mesh_manager_free()
{
  if (ib_mesh_manager.meshes != NULL) {
    for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
      ib_mesh_free(&ib_mesh_manager.meshes[mi]);
    }
    free(ib_mesh_manager.meshes);
    ib_mesh_manager.nm = 0;
  }
}

/**
 * @brief Generate the stencils for each node of a given mesh.
 */
void
generate_stencil_cache(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].stencil.n = 0;
    const double delta = (L0 / (1 << grid->maxdepth));
#if dimension == 2
    for (int mi = -2; mi <= 2; mi++) {
      for (int mj = -2; mj <= 2; mj++) {
        Point point =
          locate(POS_PBC_X(mesh->nodes[ni].lagpos.x + mi * delta),
                 POS_PBC_Y(mesh->nodes[ni].lagpos.y + mj * delta));
        if (point.level >= 0 && point.level != grid->maxdepth) {
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        }

        cache_append(&(mesh->nodes[ni].stencil), point, 0);
      }
    }
#else // dimension == 3
    for (int mi = -2; mi <= 2; mi++) {
      for (int mj = -2; mj <= 2; mj++) {
        for (int mk = -2; mk <= 2; mk++) {
          Point point =
            locate(POS_PBC_X(mesh->nodes[ni].lagpos.x + mi * delta),
                   POS_PBC_Y(mesh->nodes[ni].lagpos.y + mj * delta),
                   POS_PBC_Y(mesh->nodes[ni].lagpos.z + mk * delta));
          if (point.level >= 0 && point.level != grid->maxdepth) {
            fprintf(stderr,
                    "Warning: Lagrangian stencil not fully resolved.\n");
          }
          cache_append(&(mesh->nodes[ni].stencil), point, 0);
        }
      }
    }
#endif
  }
}

/**
 * @brief Generate stencils for every node of every mesh
 */
void
generate_stencil_caches()
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    generate_stencil_cache(&IBMESH(mi));
  }
}

vector forcing[];

/**
 * @brief Spread the IBM force at the Lagrangian nodes back
 * onto the Eulerian vector field
 *
 * @param forcing Eulerian vector field to hold the result in
 * @param mesh The mesh for to gather forcing from
 */
void
spread_eulerian_forcing(vector forcing, IBMesh* m)
{
  for (int ni = 0; ni < m->nn; ni++) {
    foreach_cache(m->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        coord dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, m->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, m->nodes[ni].lagpos.y);

        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / (sq(4 * Delta));
          foreach_dimension()
          {
            forcing.x[] += weight * m->nodes[ni].force.x;
          }
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, m->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, m->nodes[ni].lagpos.y);
        dist.z = GENERAL_1DIST(z, m->nodes[ni].lagpos.z);

        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          (cube(4 * Delta));
          foreach_dimension()
          {
            forcing.x[] += weight * m->nodes[ni].force.x;
          }
        }
#endif
      }
    }
  }
}

void
interpolate_eulerian_velocities(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_dimension()
    {
      mesh->nodes[ni].eulvel.x = 0;
    }
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        coord dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / 16.;
          foreach_dimension()
          {
            mesh->nodes[ni].eulvel.x += weight * u.x[];
          }
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].lagpos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) / 64.;
          foreach_dimension()
          {
            mesh->nodes[ni].eulvel.x += weight * u.x[];
          }
        }
#endif
      }
    }
  }
}

scalar stencils[];

/**
 * @brief
 *
 * @param mesh The immersed boundary mesh
 */
void tag_stencil(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        coord dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          stencils[] = sq(dist.x + dist.y) / sq(2. * Delta) * (2. + noise());
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].lagpos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          stencils[] =
            sq(dist.x + dist.y + dist.z) / cube(2. * Delta) * (2. + noise());
        }
#endif
      }
    }
  }
}

/**
 * @brief
 *
 */

void
tag_stencils()
{
  foreach () {
    stencils[] = 0.;
  }
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {
      tag_stencil(&ib_mesh_manager.meshes[mi]);
    }
  }
}

/**
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 */

void
advance_lagrangian_mesh(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 2
    mesh->nodes[ni].lagvel.x = -1.;
    mesh->nodes[ni].lagvel.y = 0.0;
#else // dimension == 3
    mesh->nodes[ni].lagvel.x = -1.;
    mesh->nodes[ni].lagvel.y = 0.0;
    mesh->nodes[ni].lagvel.z = 0.0;
#endif
    foreach_dimension()
    {
      mesh->nodes[ni].lagpos.x =
        mesh->nodes[ni].lagpos.x + dt * mesh->nodes[ni].lagvel.x;
    }
  }
}

/**
 * @brief Update the residual between the interpolated Eulerian velocity and
 * the Lagrangian velocity at the nodes
 */
// void
// compute_kinematic_constraint_residual(IBMesh* mesh)
// {
//   for (int ni = 0; ni < mesh->nn; ni++) {
//     mesh->nodes[ni].residual = 0;
//     foreach_dimension()
//     {
//       mesh->nodes[ni].residual += sq(mesh->nodes[ni].lagvel.x - mesh->nodes[ni].eulvel.x);
//     }
//     mesh->nodes[ni].residual = sqrt(mesh->nodes[ni].residual);
//   }
// }

void compute_kinematic_constraint_force(IBMesh *mesh) {
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_dimension() {
      mesh->nodes[ni].force.x =  (mesh->nodes[ni].lagvel.x - mesh->nodes[ni].eulvel.x);
    }
  }
}

/**
 * @brief
 * 
 * See [centered](https://basilisk.fr/src/navier-stokes/centered.h#time-integration).
 */

event
defaults(i = 0)
{
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
    {
      a.x[] = 0.;
    }
  }
}

/**
 * @brief
 *
 * See
 * [centered](https://basilisk.fr/src/navier-stokes/centered.h#time-integration).
 * Our event here will overload and run before the previously defined projection
 * event, giving us access to the predicted \f(\vec{u}^*\f).
 */

event
tracer_advection(i++)
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {

      // Advance the mesh in time
      advance_lagrangian_mesh(&ib_mesh_manager.meshes[mi]);

      // Generate the stencil cache for the nodes and the current grid
      generate_stencil_cache(&ib_mesh_manager.meshes[mi]);

      // Reset lagrangian force
      for (int ni = 0; ni < ib_mesh_manager.nm; ni++) {
        foreach_dimension()
        {
          IBMESH(mi).nodes[ni].force.x = 0.;
        }
      }
    }
  }
}

/**
 * @brief
 *
 * See
 * [centered](https://basilisk.fr/src/navier-stokes/centered.h#acceleration-term)
 */
event
acceleration(i++)
{
  // Set the previous forcing to zero
  face vector ae = a;
  foreach () {
    if (cm[] > 1e-20) {
      foreach_dimension()
      {
        forcing.x[] = 0.;
      }
    }
  }

  // Interpolate velocities
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {
      interpolate_eulerian_velocities(&ib_mesh_manager.meshes[mi]);
    }
  }

  // Compute the forcing from kinematic condition (e.g. no-slip)
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {
      compute_kinematic_constraint_force(&ib_mesh_manager.meshes[mi]);
    }
  }

  //
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {
      spread_eulerian_forcing(forcing, &ib_mesh_manager.meshes[mi]);
    }
  }

  // Contribute forcing to our acceleration face field
  foreach_face()
  {
    if (fm.x[] > 1e-20) {
      ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
    }
  }
}

/**
 * @brief
 *
 */


/**
 * @brief
 */
event cleanup(t = end) {
  ib_mesh_manager_free();
}

/**
 * Here we set up the main program
 */

#include "embed.h"
#include "tracer.h"


scalar f[];
scalar * tracers = {f};
double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];

double U0 = 1.;

int main()
{
  L0 = 8;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
  
  run(); 
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);


event init (t = 0)
{  
  foreach()
    u.x[] = cs[] ? U0 : 0.;

  ib_mesh_manager_init(1);

  ib_mesh_init(&ib_mesh_manager.meshes[0], 2);

  ib_mesh_manager.meshes->nodes[0].lagpos.x =  6.0;
  ib_mesh_manager.meshes->nodes[0].lagpos.y =  0.1;
  ib_mesh_manager.meshes->nodes[1].lagpos.x =  6.0;
  ib_mesh_manager.meshes->nodes[1].lagpos.y = -0.1;
}

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

#include "library/output_htg.h"

event movies (i += 4; t <= 10.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
  output_hdf_htg({ omega, p, stencils }, { u }, "ex2");
}


event adapt (i++) {
  tag_stencils();
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
  generate_stencil_caches();
}