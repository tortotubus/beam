#pragma once

#define BGHOSTS 2
#define DEBUG 1

#include "interface/ibm/immersedboundary.h"
#include "library/navier-stokes/centered.h"

/**
 * Distance functions
 */

#include "library/ibm/kernels.h"

/**
 * IBNode
 */
#if dimension == 1
#define STENCIL_SIZE 5
#elif dimension == 2
#define STENCIL_SIZE 25
#else // dimension == 3
#define STENCIL_SIZE 125
#endif

/**
 * @brief Initialize an IBNode struct
 */
typedef struct
{
  vertex_t lagpos;
  vertex_t force;
  vertex_t eulvel;
  vertex_t lagvel;

  // vertex_t pos;
  // vertex_t velocity;
  // vertex_t force;

  Cache stencil;
} IBNode;

void
ibnode_init_stencil(IBNode* n)
{
  n->stencil.n = 0;
  n->stencil.nm = STENCIL_SIZE;
  n->stencil.p = (Index*)malloc(STENCIL_SIZE * sizeof(Index));
}

void
ibnode_free_stencil(IBNode* n)
{
  free(n->stencil.p);
}

void
ibnode_init(IBNode* n)
{
  ibnode_init_stencil(n);
}

void
ibnode_free(IBNode* n)
{
  ibnode_free_stencil(n);
}

/**
 * IBMesh
 */

typedef enum
{
  INVALID = -1,
  NONE = 0,
  VELOCITY_COUPLED_IBM = 1,
  FORCE_COUPLED_IBM = 2
} IBMeshModelType;

typedef union
{
  ib_velocity_structure_model_t velocity_coupled;
  ib_force_structure_model_t force_coupled;
} IBMeshModel;

typedef struct
{
  IBMeshModelType model_type;
  IBMeshModel model;
  IBNode* nodes;
  int nn;
  bool isactive;
} IBMesh;

void
ib_mesh_free(IBMesh* mesh)
{
  if (mesh->nodes) {
    for (int i = 0; i < mesh->nn; i++) {
      ibnode_free(&mesh->nodes[i]);
    }
    free(mesh->nodes);
  }
}

/**
 * @brief Initialize a mesh without nodes.
 *
 * @param mesh Mesh to initialize
 */
void
ib_mesh_init(IBMesh* mesh)
{
  mesh->model_type = NONE;
  mesh->nodes = NULL;
  mesh->nn = 0;
  mesh->isactive = false;
}

/**
 * @brief Initialize the IBMesh with a number of IBNodes
 *
 * @param mesh The IBMesh to initialize
 * @param nn The number of nodes to initialize
 */
void
ib_mesh_init_nn(IBMesh* mesh, int nn)
{
  // mesh->model_type = NONE;
  mesh->nodes = (IBNode*)calloc(nn, sizeof(IBNode));
  mesh->nn = nn;
  mesh->isactive = true;

  for (int i = 0; i < nn; i++) {
    ibnode_init(&mesh->nodes[i]);
  }
}

/**
 * @brief Set the (velocity-coupled) constitutive model for a given IBMesh
 *
 * @param mesh Pointer (reference) to the IBMesh we want to set the Lagrangian
 * model for
 * @param model The C++ velocity-coupled to use for this mesh
 */
void
ib_mesh_set_velocity_model(IBMesh* mesh, ib_velocity_structure_model_t model)
{
  // Free the previous mesh, if any
  ib_mesh_free(mesh);

  // Set the model pointer
  mesh->model.velocity_coupled = model;

  // Set the model type
  mesh->model_type = VELOCITY_COUPLED_IBM;

  // Get the number of nodes from the C++ model
  int nn = ib_velocity_structure_model_get_number_of_nodes(model);

  // Initialize the basilisk C mesh
  ib_mesh_init_nn(mesh, nn);

  // Get the current-time mesh from the C++ model
  ib_structure_mesh_t cpp_mesh = ib_velocity_structure_model_get_current(model);

  // Copy the points to the basilisk C mesh
  for (int i = 0; i < nn; i++) {
    mesh->nodes[i].lagpos = cpp_mesh.position[i];
  }

  // Free the C++ mesh
  ib_structure_mesh_free(&cpp_mesh);
}

void
ib_mesh_set_force_model(IBMesh* mesh, ib_force_structure_model_t model)
{
  // Free the previous mesh, if any
  ib_mesh_free(mesh);

  // Set the model pointer
  mesh->model.force_coupled = model;

  // Set the model type
  mesh->model_type = FORCE_COUPLED_IBM;

  // Get the number of nodes from the C++ model
  int nn = ib_force_structure_model_get_number_of_nodes(model);

  // Initialize the basilisk C mesh
  ib_mesh_init_nn(mesh, nn);

  // Get the current-time mesh from the C++ model
  ib_structure_mesh_t cpp_mesh = ib_force_structure_model_get_current(model);

  // Copy the points to the basilisk C mesh
  for (int i = 0; i < nn; i++) {
    mesh->nodes[i].lagpos.x = cpp_mesh.position[i].x;
    mesh->nodes[i].lagpos.y = cpp_mesh.position[i].y;
    mesh->nodes[i].lagpos.z = cpp_mesh.position[i].z;
  }

  // Free the C++ mesh
  ib_structure_mesh_free(&cpp_mesh);
}

/**
 * @brief Helper function to compute the centroid of a given IBMesh
 *
 * @param mesh IBMesh to compute the centroid of
 */
vertex_t
ib_mesh_get_centroid(IBMesh* mesh)
{
  vertex_t centroid = { .x = 0, .y = 0, .z = 0 };
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 2
    centroid.x += mesh->nodes[ni].lagpos.x;
    centroid.y += mesh->nodes[ni].lagpos.y;
#else // dimension == 3
    centroid.x += mesh->nodes[ni].lagpos.x;
    centroid.y += mesh->nodes[ni].lagpos.y;
    centroid.z += mesh->nodes[ni].lagpos.z;
#endif
  }

  if (centroid.x != 0)
    centroid.x /= (double)mesh->nn;
  if (centroid.y != 0)
    centroid.y /= (double)mesh->nn;
  if (centroid.z != 0)
    centroid.z /= (double)mesh->nn;

  return centroid;
}

/**
 * @brief
 */
void
ib_mesh_print_pos(IBMesh* mesh)
{
  printf("{");
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 2
    if (ni == mesh->nn - 1) {
      printf("(%f,%f)\n", mesh->nodes[ni].lagpos.x, mesh->nodes[ni].lagpos.y);
    } else {
      printf("(%f,%f), ", mesh->nodes[ni].lagpos.x, mesh->nodes[ni].lagpos.y);
    }
#else // dimension == 3
    if (ni == mesh->nn - 1) {
      printf("(%f,%f,%f)\n",
             mesh->nodes[ni].lagpos.x,
             mesh->nodes[ni].lagpos.y,
             mesh->nodes[ni].lagpos.z);
    } else {
      printf("(%f,%f,%f), ",
             mesh->nodes[ni].lagpos.x,
             mesh->nodes[ni].lagpos.y,
             mesh->nodes[ni].lagpos.z);
    }
#endif
  }
  printf("}\n");
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
  if (ib_mesh_manager.meshes != NULL) {
    for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
      ib_mesh_init(&IBMESH(mi));
    }
  }

  ib_mesh_manager.nm = nm;
  ib_mesh_manager.meshes = (IBMesh*)calloc(nm, sizeof(IBMesh));

  for (int mi = 0; mi < nm; mi++) {
    ib_mesh_init(&IBMESH(mi));
  }
}

/**
 * @brief Free all members in the immersed boundary mesh manager, including the
 * array of meshes, their nodes, and stencil nodes.
 */
void
ib_mesh_manager_free()
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    ib_mesh_free(&IBMESH(mi));
  }

  ib_mesh_manager.nm = 0;
}

/**
 * @brief Generate the stencils for each node of a given mesh.
 *
 * @param mesh
 */
void
generate_stencil_cache(IBMesh* mesh)
{
  for (size_t ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].stencil.n = 0;
    double delta = (L0 / (1 << grid->maxdepth));
#if dimension == 2
    for (int mi = -2; mi <= 2; mi++) {
      for (int mj = -2; mj <= 2; mj++) {
        Point point = locate(POS_PBC_X(mesh->nodes[ni].lagpos.x + mi * delta),
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
                   POS_PBC_Z(mesh->nodes[ni].lagpos.z + mk * delta));
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
trace void
generate_stencil_caches()
{
  for (int mi = 0; mi < ib_mesh_manager.nm; ++mi) {
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
trace void
spread_eulerian_forcing(vector forcing, IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / (sq(4 * Delta));

          foreach_dimension()
          {
            forcing.x[] += weight * mesh->nodes[ni].force.x;
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
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          (cube(4 * Delta));

          foreach_dimension()
          {
            forcing.x[] += weight * mesh->nodes[ni].force.x;
          }
        }
#endif
      }
    }
  }
}
/**
 * @brief Interpolate Eulerian velocities from the Eulerian mesh to the
 * Lagrangian nodes
 *
 * @param mesh The mesh for to gather forcing from
 */
trace void
interpolate_eulerian_velocities(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_dimension()
    {
      mesh->nodes[ni].eulvel.x = 0.;
    }
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // Ignore
#elif dimension == 2
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
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].pos.z);
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
trace void
tag_stencil(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
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
trace void
tag_stencils()
{
  foreach () {
    stencils[] = 0.;
  }

  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (IBMESH(mi).isactive) {
      tag_stencil(&IBMESH(mi));
    }
  }
}

/**
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 * @param mesh
 */
trace void
advect_lagrangian_mesh(IBMesh* mesh)
{
  assert(mesh->model_type == 1); // VELOCITY_COUPLED_IBM
#if 1                            // First order
  vertex_t* velocity = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    velocity[ni].x = mesh->nodes[ni].eulvel.x;
    velocity[ni].y = mesh->nodes[ni].eulvel.y;
    velocity[ni].z = mesh->nodes[ni].eulvel.z;
    // printf("velocity[%d].x = %f\n", ni, velocity[ni].x);
    // printf("mesh->nodes[%d].velocity.x = %f\n", ni,
    // mesh->nodes[ni].velocity.x);
  }

  ib_structure_mesh_t vel_str_mesh = ib_velocity_structure_model_get_next(
    mesh->model.velocity_coupled, velocity, mesh->nn, dt);

  free(velocity);

  for (int ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].lagpos.x = vel_str_mesh.position[ni].x;
    mesh->nodes[ni].lagpos.y = vel_str_mesh.position[ni].y;
    mesh->nodes[ni].lagpos.z = vel_str_mesh.position[ni].z;
  }

  ib_structure_mesh_free(&vel_str_mesh);

#else // Second order

#endif
}

/**
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 * @param mesh
 */
void
advance_lagrangian_mesh(IBMesh* mesh)
{
#if DEBUG
  printf("advance_lagrangian_mesh\n");
#endif

  assert(mesh->model_type == 2); // FORCE_COUPLED_MESH

  vertex_t* force = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    force[ni].x = -mesh->nodes[ni].force.x;
    force[ni].y = -mesh->nodes[ni].force.y;
    force[ni].z = -mesh->nodes[ni].force.z;
  }
  
  ib_structure_mesh_t cpp_mesh = ib_force_structure_model_get_next(mesh->model.force_coupled, force, mesh->nn, dt);

  free(force);

  for (int ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].lagpos.x = cpp_mesh.position[ni].x;
    mesh->nodes[ni].lagpos.y = cpp_mesh.position[ni].y;
    mesh->nodes[ni].lagpos.z = cpp_mesh.position[ni].z;
    mesh->nodes[ni].lagvel.x = cpp_mesh.velocity[ni].x;
    mesh->nodes[ni].lagvel.y = cpp_mesh.velocity[ni].y;
    mesh->nodes[ni].lagvel.z = cpp_mesh.velocity[ni].z;
  }

  ib_structure_mesh_free(&cpp_mesh);

//   for (int ni = 0; ni < mesh->nn; ni++) {
// #if dimension == 2
//     mesh->nodes[ni].lagvel.x = -1.;
//     mesh->nodes[ni].lagvel.y = 0.0;
// #else // dimension == 3
//     mesh->nodes[ni].lagvel.x = -1.;
//     mesh->nodes[ni].lagvel.y = 0.0;
//     mesh->nodes[ni].lagvel.z = 0.0;
// #endif
//     foreach_dimension()
//     {
//       mesh->nodes[ni].lagpos.x =
//         mesh->nodes[ni].lagpos.x + dt * mesh->nodes[ni].lagvel.x;
//     }
//   }
}

/**
 * @brief Update the force between the interpolated Eulerian velocity and
 * the Lagrangian velocity at the nodes
 *
 * @param mesh
 */
void
compute_kinematic_constraint_force(IBMesh* mesh)
{
  #if DEBUG
  printf("compute_kinematic_constraint_force\n");
  #endif

  assert(mesh->model_type == 2); // FORCE_COUPLED_IBM

  double force_norm_val = 0;

  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_dimension()
    {
      mesh->nodes[ni].force.x = (mesh->nodes[ni].lagvel.x - mesh->nodes[ni].eulvel.x);

      force_norm_val += mesh->nodes[ni].force.x;
    }
  }


  // printf("||f_ib|| = %f\n", force_norm_val);
}

/**
 * @brief
 *
 * See
 * [centered](https://basilisk.fr/src/navier-stokes/centered.h#time-integration).
 */
event
defaults(i = 0)
{
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

/**
 * @brief
 *
 * See
 * [centered](https://basilisk.fr/src/navier-stokes/centered.h#time-integration).
 * Our event here will overload and run before the previously defined
 * tracer_advection event.
 */
event
tracer_advection(i++)
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (IBMESH(mi).isactive) {
      // Interpolate eulerain velocities
      interpolate_eulerian_velocities(&IBMESH(mi));

      // Advance the mesh in time
      switch (IBMESH(mi).model_type) {
        case VELOCITY_COUPLED_IBM: // VELOCITY_COUPLED_IBM:
          advect_lagrangian_mesh(&IBMESH(mi));
          break;

        case FORCE_COUPLED_IBM: // FORCE_COUPLED_IBM:
          advance_lagrangian_mesh(&IBMESH(mi));
          break;

        case 0: // NONE:
          break;

        default:
          assert(1 == 2);
          break;
      }

      // Generate the stencil cache for the nodes and the current grid
      generate_stencil_cache(&IBMESH(mi));

      // Reset lagrangian force
      for (int j = 0; j < IBMESH(mi).nn; j++) {
        foreach_dimension()
        {
          // IBMESH(mi).nodes[j].force.x = 0.;
        }
      }
    }
  }
}

/**
 * @brief In the acceleration event, we transfer the Lagrangian forces to the
 * fluid using a regularized Dirac function. The acceleration is stored on the
 * cell faces, and will be fed as a source term to the Navier-Stokes solver.
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
    if (cm[] > 1.e-20) {
      foreach_dimension()
      {
        forcing.x[] = 0.;
      }
    }
  }

  // Compute the internal structural forces
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (IBMESH(mi).isactive) {

      // switch (IBMESH(mi).model_type) {
      //   case FORCE_COUPLED_IBM: {
      //     compute_kinematic_constraint_force(&IBMESH(mi));
      //     break;
      //   }
      // }

      if (IBMESH(mi).model_type == VELOCITY_COUPLED_IBM) {
        vertex_t* velocity = calloc((size_t)IBMESH(mi).nn, sizeof(vertex_t));
        for (int ni = 0; ni < IBMESH(mi).nn; ni++) {
          velocity[ni].x = IBMESH(mi).nodes[ni].eulvel.x;
          velocity[ni].y = IBMESH(mi).nodes[ni].eulvel.y;
          velocity[ni].z = IBMESH(mi).nodes[ni].eulvel.z;
        }

        ib_structure_mesh_t mesh_midpoint =
          ib_velocity_structure_model_get_midpoint(
            IBMESH(mi).model.velocity_coupled, velocity, IBMESH(mi).nn, dt);

        free(velocity);

        for (int ni = 0; ni < IBMESH(mi).nn; ni++) {
          IBMESH(mi).nodes[ni].force.x = mesh_midpoint.forces[ni].x;
          IBMESH(mi).nodes[ni].force.y = mesh_midpoint.forces[ni].y;
          IBMESH(mi).nodes[ni].force.z = mesh_midpoint.forces[ni].z;
        }

        ib_structure_mesh_free(&mesh_midpoint);
      } else {
        compute_kinematic_constraint_force(&IBMESH(mi));
      }
    }
  }

  // Spread Eulerian forcing
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (ib_mesh_manager.meshes[mi].isactive) {
      spread_eulerian_forcing(forcing, &ib_mesh_manager.meshes[mi]);
    }
  }

  foreach_face()
  {
    if (fm.x[] > 1.e-20)
      ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
  }
}

/** At the end of the simulation, we free the allocated memory.*/
event
cleanup(t = end)
{
  ib_mesh_manager_free();
}