#pragma once

#define BGHOSTS 2

#ifndef DEBUG
#define DEBUG 0
#endif

#include "interface/ibm/immersedboundary.h"
#include "library/navier-stokes/centered.h"

/* ====================================================================================================================
 * Distance functions
 * ====================================================================================================================
 */

#include "library/ibm/kernels.h"

/* ====================================================================================================================
 * IBNode
 * ====================================================================================================================
 */

#if dimension == 1
#define STENCIL_SIZE 5
#elif dimension == 2
#define STENCIL_SIZE 25
#else // dimension == 3
#define STENCIL_SIZE 125
#endif

#define RHO_F 1

/**
 * @brief Initialize an IBNode struct
 */
typedef struct
{
  vertex_t lagpos; /**<  lagrangian nodal position */
  vertex_t force;  /**< lagrange multiplier as a force */
  vertex_t lagvel; /**< lagrangian nodal velocity */
  vertex_t
    eulvel; /**< eulerian fluid velocity interpolated at the lagrangian node */

  vertex_t rhs; /**< b_i = lagvel_i - (Ju*)_i (constraint mismatch)*/
  vertex_t res; /**< r_i (residual) */
  vertex_t w;   /**< w_i (search direction) */
  vertex_t Ay;  /**< y_i = (A w)_i */

  double weight;    /**< quadtrature weight */
  vertex_t gravity; /**< gravity */

  Cache stencil; /**< stencil cache */
} IBNode;

/**
 * @brief Constructor for the stencil member
 *
 * @memberof IBNode
 */
void
ibnode_init_stencil(IBNode* node)
{
  node->stencil.n = 0;
  node->stencil.nm = STENCIL_SIZE;
  node->stencil.p = (Index*)malloc(STENCIL_SIZE * sizeof(Index));
}

/**
 * @brief Desctructor for the stencil member
 * @memberof IBNode
 */
void
ibnode_free_stencil(IBNode* node)
{
  free(node->stencil.p);
}

/**
 * @brief Constructor for the IBNode struct
 * @memberof IBNode
 */
void
ibnode_init(IBNode* node)
{
  foreach_dimension()
  {
    node->lagpos.x = 0;
    node->force.x = 0;
    node->lagvel.x = 0;
    node->eulvel.x = 0;
    node->rhs.x = 0;
    node->res.x = 0;
    node->w.x = 0;
    node->Ay.x = 0;
    node->gravity.x = 0;
  }

  node->weight = 1.;

  ibnode_init_stencil(node);
}

/**
 * @brief Destructor for the IBNode struct
 * @memberof IBNode
 */
void
ibnode_free(IBNode* node)
{
  ibnode_free_stencil(node);
}

/* ====================================================================================================================
 * IBMesh
 * ====================================================================================================================
 */

typedef struct
{
  ib_force_structure_model_t model;
  IBNode* nodes;
  int nn;
  bool isactive;
} IBMesh;

/**
 * @brief Destructor for the IBMesh struct
 * @memberof IBMesh
 */
void
ib_mesh_free(IBMesh* mesh)
{
  if (mesh->nodes) {
    for (int i = 0; i < mesh->nn; i++) {
      ibnode_free(&mesh->nodes[i]);
    }
    free(mesh->nodes);
    mesh->nodes = NULL;
  }
  mesh->nn = 0;
  mesh->isactive = false;
}

/**
 * @brief Default constructor for an IBMesh
 *
 * @param mesh Mesh to initialize
 *
 * @memberof IBMesh
 */
void
ib_mesh_init(IBMesh* mesh)
{
  mesh->model = NULL;
  mesh->nodes = NULL;
  mesh->nn = 0;
  mesh->isactive = false;
}

/**
 * @brief Construct the IBMesh with a number of IBNodes
 *
 * @param mesh The IBMesh to initialize
 * @param nn The number of nodes to initialize
 *
 * @memberof IBMesh
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
 * @brief Set the model for an IBMesh
 *
 * @param mesh The IBMesh
 * @param model The pointer to the force-coupled structural model
 *
 * @memberof IBMesh
 */
void
ib_mesh_model(IBMesh* mesh, ib_force_structure_model_t model)
{
#if DEBUG
  printf("ib_mesh_model\n");
#endif
  // Free the previous mesh, if any
  ib_mesh_free(mesh);

  // Set the model pointer
  mesh->model = model;

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
 *
 * @memberof IBMesh
 */
vertex_t
ib_mesh_get_centroid(IBMesh* mesh)
{
  vertex_t centroid = { .x = 0, .y = 0, .z = 0 };
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 1
    // TODO
#elif dimension == 2
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
 * @brief Prints the centroid of an IBMesh
 *
 * @param mesh The IBMesh
 *
 * @memberof IBMesh
 */
void
ib_mesh_print_pos(IBMesh* mesh)
{
  printf("{");
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 1
    // TODO
#elif dimension == 2
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
 * @brief Generate the stencils for each node of a given mesh.
 *
 * @param mesh
 *
 * @memberof IBMeshManager
 */
trace void
ib_mesh_generate_stencil_cache(IBMesh* mesh)
{
#if DEBUG
  printf("ib_mesh_generate_stencil_cache\n");
#endif
  for (size_t ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].stencil.n = 0;
    double delta = (L0 / (1 << grid->maxdepth));
#if dimension == 1
    // TODO
#elif dimension == 2
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

vector forcing[];
vector tmp_force[];
vector tmp_vel[];

/**
 * @brief Spread the IBM force at the Lagrangian nodes back
 * onto the Eulerian vector field
 *
 * @param mesh The mesh for to gather forcing from
 * @param forcing Eulerian vector field to hold the result in
 *
 * @memberof IBMesh
 */
trace void
ib_mesh_spread_eulerian_forcing(IBMesh* mesh, vector forcing)
{
#if DEBUG
  printf("ib_mesh_spread_eulerian_forcing\n");
#endif
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_cache(node->stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // TODO
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / sq(4. * Delta);

          foreach_dimension()
          {
            forcing.x[] += weight * node->force.x * node->weight;
          }
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        dist.z = GENERAL_1DIST(z, node->lagpos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          cube(4. * Delta);

          foreach_dimension()
          {
            forcing.x[] += weight * node->force.x * node->weight;
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
 *
 * @memberof IBMesh
 */
trace void
ib_mesh_interpolate_eulerian_velocities(IBMesh* mesh)
{
#if DEBUG
  printf("ib_mesh_interpolate_eulerian_velocities\n");
#endif
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_dimension()
    {
      mesh->nodes[ni].eulvel.x = 0.;
    }

    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // TODO
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / sq(4. * Delta);
          foreach_dimension()
          {
            node->eulvel.x += weight * u.x[];
          }
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        dist.z = GENERAL_1DIST(z, node->lagpos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          cube(4. * Delta);
          foreach_dimension()
          {
            node->eulvel.x += weight * u.x[];
          }
        }
#endif
      }
    }

#if _MPI
    for (int ni = 0; ni < mesh->nn; ni++) {
      IBNode* node = &mesh->nodes[ni];
      double buf[3] = { 0., 0., 0. };
#if dimension == 1
      buf[0] = node->eulvel.x;
#elif dimension == 2
      buf[0] = node->eulvel.x;
      buf[1] = node->eulvel.y;
#else // dimension == 3
      buf[0] = node->eulvel.x;
      buf[1] = node->eulvel.y;
      buf[2] = node->eulvel.z;
#endif

      MPI_Allreduce(
        MPI_IN_PLACE, buf, dimension, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#if dimension == 1
      node->eulvel.x = buf[0];
#elif dimension == 2
      node->eulvel.x = buf[0];
      node->eulvel.y = buf[1];
#else // dimension == 3
      node->eulvel.x = buf[0];
      node->eulvel.y = buf[1];
      node->eulvel.z = buf[2];
#endif
    }
#endif
  }
}

scalar stencils[];

/**
 * @brief Injects the stencils field with noise so that we may call adapt_wavlet
 * to ensure sufficient cells on the same level near our IBNodes
 *
 * @param mesh The immersed boundary mesh
 *
 * @memberof IBMesh
 */
trace void
ib_mesh_tag_stencil(IBMesh* mesh)
{
#if DEBUG
  printf("tag_stencil\n");
#endif
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_cache(node->stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // TODO
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          stencils[] = sq(dist.x + dist.y) / sq(2. * Delta) * (2. + noise());
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        dist.z = GENERAL_1DIST(z, node->lagpos.z);
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
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 * @param mesh
 *
 * @memberof IBMesh
 */
trace void
ib_mesh_advance_lagrangian_mesh(IBMesh* mesh)
{
#if DEBUG
  printf("ib_mesh_advance_lagrangian_mesh\n");
#endif

  vertex_t* force = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    force[ni].x = -node->force.x * node->weight + node->gravity.x;
    force[ni].y = -node->force.y * node->weight + node->gravity.y;
    force[ni].z = -node->force.z * node->weight + node->gravity.z;
  }

  ib_structure_mesh_t cpp_mesh =
    ib_force_structure_model_get_next(mesh->model, force, mesh->nn, dt);

  free(force);

  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    node->lagpos.x = cpp_mesh.position[ni].x;
    node->lagpos.y = cpp_mesh.position[ni].y;
    node->lagpos.z = cpp_mesh.position[ni].z;
    node->lagvel.x = cpp_mesh.velocity[ni].x;
    node->lagvel.y = cpp_mesh.velocity[ni].y;
    node->lagvel.z = cpp_mesh.velocity[ni].z;
  }

  ib_structure_mesh_free(&cpp_mesh);
}

/**
 * @brief Compute the kinematic no-slip residual
 * 
 * @memberof IBMesh
 */
trace void
ib_mesh_compute_constraint_rhs(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_dimension()
    {
      node->rhs.x = node->lagvel.x - node->eulvel.x; // b = v_lag - J u*
    }
  }
}

/**
 * @memberof IBMesh
 */
trace void
ib_mesh_matvec_Aw(IBMesh* mesh, double dt)
{

  foreach () {
    foreach_dimension()
    {
      tmp_force.x[] = 0;
      tmp_vel.x[] = 0;
    }
  }

  // 1. g = J^T w on grid
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];

    foreach_cache(node->stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;

#if dimension == 1
        // TODO
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {

          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / sq(4. * Delta);

          foreach_dimension()
          {
            tmp_force.x[] += weight * node->w.x * node->weight;
          }
        }
#else
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].lagpos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].lagpos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].lagpos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {

          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          cube(4. * Delta);

          foreach_dimension()
          {
            tmp_force.x[] += weight * node->w.x * node->weight;
          }
        }
#endif
      }
    }
  }

  // 2. c = M_f^{-1} g ~ g / (rho_f)
  foreach ()
    foreach_dimension() tmp_vel.x[] = tmp_force.x[] / RHO_F;

  // 3. d = J c, store in node->Ay temporarily, then scale by dt
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];

    // Clear the vector row
    foreach_dimension()
    {
      node->Ay.x = 0.;
    }

    foreach_cache(node->stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // TODO
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);

        if (fabs(dist.x) <= 2. * Delta && fabs(dist.y) <= 2. * Delta) {
          double weight = (1. + cos(0.5 * pi * dist.x / Delta)) *
                          (1. + cos(0.5 * pi * dist.y / Delta)) /
                          sq(4. * Delta);

          foreach_dimension()
          {
            node->Ay.x += weight * tmp_vel.x[];
          }
        }
#else
        dist.x = GENERAL_1DIST(x, node->lagpos.x);
        dist.y = GENERAL_1DIST(y, node->lagpos.y);
        dist.z = GENERAL_1DIST(z, node->lagpos.z);

        if (fabs(dist.x) <= 2. * Delta && fabs(dist.y) <= 2. * Delta &&
            fabs(dist.z) <= 2. * Delta) {
          double weight = (1. + cos(0.5 * pi * dist.x / Delta)) *
                          (1. + cos(0.5 * pi * dist.y / Delta)) *
                          (1. + cos(0.5 * pi * dist.z / Delta)) /
                          cube(4. * Delta);

          foreach_dimension()
          {
            node->Ay.x += weight * tmp_vel.x[] * node->weight;
          }
        }
#endif
      }
    }

    // Final scaling y = A w = Δt * J M_f^{-1} J^T w
    foreach_dimension()
    {
      node->Ay.x *= dt;
    }
  }
}

/**
 * @memberof IBMesh
 */
trace void
ib_mesh_solve_lambda_CG(IBMesh* mesh, double dt)
{
  const int maxiter = 500;
  const double tol = 1e-6;

  // Initial guess
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_dimension()
    {
      node->force.x = 0.;        // lambda^0 = 0
      node->res.x = node->rhs.x; // res^0 = b
      node->w.x = node->res.x;   // w^0 = res^0
    }
  }

  double rr_old = 0.;
  for (int ni = 0; ni < mesh->nn; ni++) {
    IBNode* node = &mesh->nodes[ni];
    foreach_dimension()
    {
      rr_old += sq(node->res.x) * node->weight;
    }
  }

  // printf("=======================================================\n");
  if (sqrt(rr_old) < tol) {
    return;
  }

  for (int k = 0; k < maxiter; k++) {
    // printf("k: %d : ||R|| = %f\n", k, rr_old);
    ib_mesh_matvec_Aw(mesh, dt);

    double wy = 0.;
    for (int ni = 0; ni < mesh->nn; ni++) {
      IBNode* node = &mesh->nodes[ni];
      foreach_dimension()
      {
        wy += node->w.x * node->Ay.x * node->weight;
      }
    }

    if (fabs(wy) < 1e-30) {
      break;
    }
    double alpha = rr_old / wy;
    double rr_new = 0.;

    // \f(\lambda_{k+1} = \lambda_{k} + \alpha w_{k}\f), and \f(r_{k+1} = r_{k}
    // - \alpha y_k\f).
    for (int ni = 0; ni < mesh->nn; ni++) {
      IBNode* node = &mesh->nodes[ni];
      foreach_dimension()
      {
        node->force.x += alpha * node->w.x; // lambda update
        node->res.x -= alpha * node->Ay.x;  // residual update
        rr_new += sq(node->res.x) * node->weight;
      }
    }

    // Check tolerance
    if (sqrt(rr_new) < tol) {
      break;
    }

    double beta = rr_new / rr_old;
    rr_old = rr_new;

    // w_{k+1} = r_{k+1} + \beta w_{k}
    for (int ni = 0; ni < mesh->nn; ni++) {
      IBNode* node = &mesh->nodes[ni];
      foreach_dimension()
      {
        node->w.x = node->res.x + beta * node->w.x;
      }
    }
  }
}

/* ====================================================================================================================
 * IBMeshManager
 * ====================================================================================================================
 */

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
 * @brief Free all members in the immersed boundary mesh manager, including the
 * array of meshes, their nodes, and stencil nodes.
 *
 * @memberof IBMeshManager
 */
void
ib_mesh_manager_free()
{
#if DEBUG
  printf("ib_mesh_manager_free\n");
#endif
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    ib_mesh_free(&IBMESH(mi));
  }

  free(ib_mesh_manager.meshes);
  ib_mesh_manager.meshes = NULL;
  ib_mesh_manager.nm = 0;
}

/**
 * @brief Initialize the immersed boundary mesh manager
 *
 * @param nm The number of meshes you plan to have
 *
 * @memberof IBMeshManager
 */
void
ib_mesh_manager_init(int nm)
{
#if DEBUG
  printf("ib_mesh_manager_init\n");
#endif
  if (ib_mesh_manager.meshes != NULL) {
    ib_mesh_manager_free();
  }

  ib_mesh_manager.nm = nm;
  ib_mesh_manager.meshes = (IBMesh*)calloc(nm, sizeof(IBMesh));

  for (int mi = 0; mi < nm; mi++) {
    ib_mesh_init(&IBMESH(mi));
  }
}

/**
 * @brief Generate stencils for every node of every mesh
 *
 * @memberof IBMeshManager
 */
trace void
ib_mesh_manager_generate_stencil_caches()
{
#if DEBUG
  printf("ib_mesh_manager_generate_stencil_caches\n");
#endif
  for (int mi = 0; mi < ib_mesh_manager.nm; ++mi) {
    ib_mesh_generate_stencil_cache(&IBMESH(mi));
  }
}

/**
 * @brief
 *
 */
trace void
ib_mesh_manager_tag_stencils()
{
#if DEBUG
  printf("tag_stencils\n");
#endif

  foreach () {
    stencils[] = 0.;
  }

  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    if (IBMESH(mi).isactive) {
      ib_mesh_tag_stencil(&IBMESH(mi));
    }
  }
}

/* ====================================================================================================================
 * Centered Navier-Stokes Solver Events
 * ====================================================================================================================
 */

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

event
acceleration(i++)
{
  face vector ae = a;

  // 0. Clear the forcing term on the Eulerian grid
  foreach () {
    if (cm[] > 1e-20) {
      foreach_dimension()
      {
        forcing.x[] = 0.;
      }
    }
  }

  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    IBMesh* mesh = &ib_mesh_manager.meshes[mi];
    if (mesh->isactive) {
      // a. Advance lagrangian mesh
      ib_mesh_advance_lagrangian_mesh(mesh);

      // b.
      ib_mesh_generate_stencil_cache(mesh);

      // 1. Interpolate \f(u*\f) with \f(J u*\f)
      ib_mesh_interpolate_eulerian_velocities(mesh);

      // 2. Build the right-hand side vector \f(b = v_{\rm lag} - J u^* \f)
      ib_mesh_compute_constraint_rhs(mesh);

      // 3. Solve A \lambda = b via CG/Uzawa on nodes
      ib_mesh_solve_lambda_CG(mesh, dt);

      // 4. Spread lambda to grid forcing \f(J^T \lambda\f)
      ib_mesh_spread_eulerian_forcing(mesh, forcing);

      // 5. Add forcing to the acceleration field
      foreach_face()
      {
        if (fm.x[] > 1e-20) {
          ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
        }
      }
    }
  }
}

/** At the end of the simulation, we free the allocated memory.*/
event
cleanup(t = end)
{
  ib_mesh_manager_free();
}

event
adapt(i++)
{
  ib_mesh_manager_tag_stencils();
  // adapt_wavelet({ stencils }, (double[]){ 1.e-2 }, maxlevel = 9);
  // ib_mesh_manager_generate_stencil_caches();
}