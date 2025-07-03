#define BGHOSTS 2

#ifndef DEBUG_PRINT_CALL_OR_EVENT
  #define DEBUG_PRINT_CALL_OR_EVENT 0
#endif

#include "interface/immersedboundary.h"
#include "library/navier-stokes/centered.h"

/**
 * Distance functions
 */

#define ACROSS_PERIODIC(a, b) (fabs(a - b) > L0 / 2.)
#define PERIODIC_1DIST(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DIST(a, b) : a - b)
#define PERIODIC_1DAVG(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DAVG(a, b) : a + b)
#define GENERAL_SQNORM(a, b) (sq(GENERAL_1DIST(a.x, b.x)) + sq(GENERAL_1DIST(a.y, b.y)) + sq(GENERAL_1DIST(a.z, b.z)))

#define POS_PBC_X(X) ((u.x.boundary[left]  != periodic_bc) ? (X) : (((X - (X0 + L0 / 2)) > L0 / 2.) ? (X) - L0 : (X)))
#define POS_PBC_Y(Y) ((u.x.boundary[top]   != periodic_bc) ? (Y) : (((Y - (Y0 + L0 / 2)) > L0 / 2.) ? (Y) - L0 : (Y)))
#define POS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : (((Z - (Z0 + L0 / 2)) > L0 / 2.) ? (Z) - L0 : (Z)))

/**
 * DFNode
 */
#if dimension == 1
#define STENCIL_SIZE 5
#elif dimension == 2
#define STENCIL_SIZE 25
#else // dimension == 3
#define STENCIL_SIZE 125
#endif

typedef struct
{
  vertex_t pos;
  vertex_t velocity;
  vertex_t force;
  Cache stencil;
} DFNode;

void
ibnode_init_stencil(DFNode* n)
{
  n->stencil.n = STENCIL_SIZE;
  n->stencil.nm = STENCIL_SIZE;
  n->stencil.p = (Index*)malloc(STENCIL_SIZE * sizeof(Index));
}

void
ibnode_free_stencil(DFNode* n)
{
  free(n->stencil.p);
}

void
ibnode_init(DFNode* n)
{
  ibnode_init_stencil(n);
}

void
ibnode_free(DFNode* n)
{
  ibnode_free_stencil(n);
}

/**
 * DFMesh
 */

typedef struct
{
  df_structure_model_t model;
  DFNode* nodes;
  int nn;
  bool isactive;
} DFMesh;

void
df_mesh_free(DFMesh* mesh)
{
  if (mesh->nodes) {
    for (int i = 0; i < mesh->nn; i++) {
      ibnode_free(&mesh->nodes[i]);
    }
    free(mesh->nodes);
  }
}

void
df_mesh_init(DFMesh* mesh, int nn)
{
  mesh->nodes = (DFNode*)calloc(nn, sizeof(DFNode));
  mesh->nn = nn;
  mesh->isactive = true;
  for (int i = 0; i < nn; i++) {
    ibnode_init(&mesh->nodes[i]);
  }
}

void
df_mesh_set_model(DFMesh* mesh, df_structure_model_t model)
{
  df_mesh_free(mesh);
  mesh->model = model;
  int nn = df_structure_model_get_number_of_nodes(model);
  df_mesh_init(mesh, nn);
  df_structure_mesh_t cpp_mesh = df_structure_model_get_current(model);
  for (int i = 0; i < nn; i++) {
    mesh->nodes[i].pos = cpp_mesh.points[i];
  }
  df_structure_mesh_free(&cpp_mesh);
}

vertex_t
df_mesh_get_centroid(DFMesh* mesh)
{
  vertex_t centroid = { .x = 0, .y = 0, .z = 0 };
  for (int i = 0; i < mesh->nn; i++) {
#if dimension == 2
    centroid.x += mesh->nodes[i].pos.x;
    centroid.y += mesh->nodes[i].pos.y;
#else // dimension == 3
    centroid.x += mesh->nodes[i].pos.x;
    centroid.y += mesh->nodes[i].pos.y;
    centroid.z += mesh->nodes[i].pos.z;
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

void
df_mesh_print_pos(DFMesh* mesh)
{
  printf("{");
  for (int ni = 0; ni < mesh->nn; ni++) {
#if dimension == 2
    if (ni == mesh->nn - 1) {
      printf("(%f,%f)", mesh->nodes[ni].pos.x, mesh->nodes[ni].pos.y);
    } else {
      printf("(%f,%f), ", mesh->nodes[ni].pos.x, mesh->nodes[ni].pos.y);
    }
#else // dimension == 3
    if (ni == mesh->nn - 1) {
      printf("(%f,%f,%f)",
             mesh->nodes[ni].pos.x,
             mesh->nodes[ni].pos.y,
             mesh->nodes[ni].pos.z);
    } else {
      printf("(%f,%f,%f), ",
             mesh->nodes[ni].pos.x,
             mesh->nodes[ni].pos.y,
             mesh->nodes[ni].pos.z);
    }
#endif
  }
  printf("}\n");
}

/**
 * IBMeshManager
 */

typedef struct
{
  DFMesh* meshes;
  int nm;
} IBMeshManager;

#define DFMesh(i) (df_mesh_manager->meshes[i])

IBMeshManager* df_mesh_manager = NULL;

void
df_mesh_manager_init(int nm)
{
  if (df_mesh_manager != NULL)
    return;
  df_mesh_manager = (IBMeshManager*)malloc(sizeof(IBMeshManager));
  df_mesh_manager->nm = nm;
  df_mesh_manager->meshes = (DFMesh*)calloc(nm, sizeof(DFMesh));
}

void
df_mesh_manager_free()
{
  if (df_mesh_manager == NULL)
    return;
  for (int i = 0; i < df_mesh_manager->nm; i++) {
    df_mesh_free(&DFMesh(i));
  }
  free(df_mesh_manager);
  df_mesh_manager = NULL;
}

void
generate_stencil(DFMesh* mesh)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: generate_stencil\n");
#endif

  for (size_t ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].stencil.n = 0;
    double delta = (L0 / (1 << grid->maxdepth));
    for (int mi = -2; mi <= 2; mi++) {
      for (int mj = -2; mj <= 2; mj++) {
#if dimension == 2
        Point point = locate(POS_PBC_X(mesh->nodes[ni].pos.x + mi * delta),
                             POS_PBC_Y(mesh->nodes[ni].pos.y + mj * delta));
#else // dimension == 3
        for (int mk = -2; mk <= 2; mk++) {
          Point point = locate(POS_PBC_X(mesh->nodes[ni].pos.x + mi * delta),
                               POS_PBC_Y(mesh->nodes[ni].pos.y + mj * delta),
                               POS_PBC_Z(mesh->nodes[ni].pos.z + mk * delta));
        }
#endif
        if (point.level >= 0 && point.level != grid->maxdepth)
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        cache_append(&(mesh->nodes[ni].stencil), point, 0);
      }
    }
  }
}

trace void
generate_stencils()
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: generate_stencils\n");
#endif
  for (int i = 0; i < df_mesh_manager->nm; ++i) {
    generate_stencil(&DFMesh(i));
  }
}

vector forcing[];

trace void
lag2eul(vector forcing, DFMesh* mesh)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: lag2eul\n");
#endif

  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / (sq(4 * Delta));

          foreach_dimension()
          {
            forcing.x[] += weight * mesh->nodes[ni].force.x;
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

trace void
eul2lag(DFMesh* mesh)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: eul2lag \n");
#endif

  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_dimension() mesh->nodes[ni].velocity.x = 0.;
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 1
        // Ignore
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / 16.;
          foreach_dimension()
          {
            mesh->nodes[ni].velocity.x += weight * u.x[];
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
            mesh->nodes[ni].velocity.x += weight * u.x[];
            printf("%f ", weight);
          }
        }
#endif
      }
    }
    // printf("\n");
  }
}

scalar stencils[];

trace void
tag_stencil(DFMesh* mesh)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: tag_stencil \n");
#endif

  for (int ni = 0; ni < mesh->nn; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        vertex_t dist;
#if dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          stencils[] = sq(dist.x + dist.y) / sq(2. * Delta) * (2. + noise());
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].pos.z);
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

trace void
tag_stencils()
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: tag_stencils \n");
#endif
  foreach () {
    stencils[] = 0.;
  }
  for (int i = 0; i < df_mesh_manager->nm; i++) {
    if (DFMesh(i).isactive) {
      tag_stencil(&DFMesh(i));
    }
  }
}

trace void
advect_mesh(DFMesh* mesh)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Function: advect_mesh\n");
#endif

  eul2lag(mesh);
#if 1 // First order
  vertex_t *velocity = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    velocity[ni].x = mesh->nodes[ni].velocity.x;
    velocity[ni].y = mesh->nodes[ni].velocity.y;
    velocity[ni].z = mesh->nodes[ni].velocity.z;
    // printf("velocity[%d].x = %f\n", ni, velocity[ni].x);
    // printf("mesh->nodes[%d].velocity.x = %f\n", ni, mesh->nodes[ni].velocity.x);
  }

  df_structure_mesh_t smesh = df_structure_model_get_next(mesh->model, velocity, mesh->nn, dt);

  free(velocity);  

  for (int ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].pos.x = smesh.points[ni].x;
    mesh->nodes[ni].pos.y = smesh.points[ni].y;
    mesh->nodes[ni].pos.z = smesh.points[ni].z;
  }

  df_structure_mesh_free(&smesh);

#else // Second order

  vertex_t* velocity = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    velocity[ni].x = mesh->nodes[ni].velocity.x;
    velocity[ni].y = mesh->nodes[ni].velocity.y;
    velocity[ni].z = mesh->nodes[ni].velocity.z;
  }

  df_structure_mesh_t smesh =
    df_structure_model_get_midpoint(mesh->model, velocity, mesh->nn, dt);

  free(velocity);

  DFMesh buffer_mesh;
  df_mesh_init(&buffer_mesh, mesh->nn);
  buffer_mesh.model = mesh->model;

  // buffer_mesh.isactive = true;
  // buffer_mesh.nn = mesh->nn;
  // buffer_mesh.nodes = calloc(mesh->nn, sizeof(DFNode));

  for (int ni = 0; ni < mesh->nn; ni++) {
    buffer_mesh.nodes[ni].pos.x = smesh.points[ni].x;
    buffer_mesh.nodes[ni].pos.y = smesh.points[ni].y;
    buffer_mesh.nodes[ni].pos.z = smesh.points[ni].z;

    buffer_mesh.nodes[ni].stencil.n = STENCIL_SIZE;
    buffer_mesh.nodes[ni].stencil.nm = STENCIL_SIZE;
    buffer_mesh.nodes[ni].stencil.p = malloc(STENCIL_SIZE * sizeof(Index));
  }

  generate_stencil(&buffer_mesh);
  eul2lag(&buffer_mesh);

  velocity = calloc(mesh->nn, sizeof(vertex_t));

  for (int ni = 0; ni < mesh->nn; ni++) {
    velocity[ni].x = buffer_mesh.nodes[ni].velocity.x;
    velocity[ni].y = buffer_mesh.nodes[ni].velocity.y;
    velocity[ni].z = buffer_mesh.nodes[ni].velocity.z;
  }

  df_structure_mesh_free(&smesh);
  smesh = df_structure_model_get_next(buffer_mesh.model, velocity, buffer_mesh.nn, dt);
  free(velocity);

  for (int ni = 0; ni < mesh->nn; ni++) {
    mesh->nodes[ni].pos.x = smesh.points[ni].x;
    mesh->nodes[ni].pos.y = smesh.points[ni].y;
    mesh->nodes[ni].pos.z = smesh.points[ni].z;
  }

  df_mesh_free(&buffer_mesh);
#endif

  generate_stencil(mesh);
}

/*
Events
*/

event
defaults(i = 0)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Event: defaults\n");
#endif

  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

event
tracer_advection(i++)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Event: tracer_advection\n");
#endif

  // printf("DFMesh(%d).isactive == %s", 0, DFMesh(0).isactive ? "true" :
  // "false");

  for (int i = 0; i < df_mesh_manager->nm; i++) {
    if (DFMesh(i).isactive) {
      advect_mesh(&DFMesh(i));
      for (int j = 0; j < DFMesh(i).nn; j++) {
        foreach_dimension()
        {
          DFMesh(i).nodes[j].force.x = 0.;
        }
      }
    }
  }
}

/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */

event
acceleration(i++)
{
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Event: acceleration\n");
#endif

  face vector ae = a;
  foreach () {
    if (cm[] > 1.e-20)
      foreach_dimension() forcing.x[] = 0.;
  }
  for (int mi = 0; mi < df_mesh_manager->nm; mi++) {
    if (DFMesh(mi).isactive) {

      // // NEW

      vertex_t* velocity = calloc((size_t)DFMesh(mi).nn, sizeof(vertex_t));
      for (int ni = 0; ni < DFMesh(mi).nn; ni++) {
        velocity[ni].x = DFMesh(mi).nodes[ni].velocity.x;
        velocity[ni].y = DFMesh(mi).nodes[ni].velocity.y;
        velocity[ni].z = DFMesh(mi).nodes[ni].velocity.z;
      }

      df_structure_mesh_t mesh_midpoint = df_structure_model_get_midpoint(
        DFMesh(mi).model, velocity, DFMesh(mi).nn, dt);
      free(velocity);

      for (int ni = 0; ni < DFMesh(mi).nn; ni++) {
        DFMesh(mi).nodes[ni].force.x = mesh_midpoint.forces[ni].x;
        DFMesh(mi).nodes[ni].force.y = mesh_midpoint.forces[ni].y;
        DFMesh(mi).nodes[ni].force.z = mesh_midpoint.forces[ni].z;
      }

      df_structure_mesh_free(&mesh_midpoint);
      // // END NEW

      lag2eul(forcing, &DFMesh(mi));
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
#if DEBUG_PRINT_CALL_OR_EVENT
  printf("Event: cleanup");
#endif

  // free_all_caps(&allCaps);
  df_mesh_manager_free();
}