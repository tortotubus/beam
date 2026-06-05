#pragma once

#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBKernels.h"
#include "elff/c/models/capsule/IBCapsule.h"

vector ib_viscosity_grid_gradient[];
scalar ib_viscosity_divergence[];

void ib_capsule_viscosity_set_wall_boundary_conditions (scalar indicator) {
  if (u.x.boundary[left] != periodic_bc) {
    indicator[left] = dirichlet (0.);
    indicator[right] = dirichlet (0.);
    ib_viscosity_divergence[left] = dirichlet (0.);
    ib_viscosity_divergence[right] = dirichlet (0.);
  }
#if dimension >= 2
  if (u.y.boundary[top] != periodic_bc) {
    indicator[top] = dirichlet (0.);
    indicator[bottom] = dirichlet (0.);
    ib_viscosity_divergence[top] = dirichlet (0.);
    ib_viscosity_divergence[bottom] = dirichlet (0.);
  }
#endif
#if dimension >= 3
  if (u.z.boundary[front] != periodic_bc) {
    indicator[front] = dirichlet (0.);
    indicator[back] = dirichlet (0.);
    ib_viscosity_divergence[front] = dirichlet (0.);
    ib_viscosity_divergence[back] = dirichlet (0.);
  }
#endif
}

static coord ib_capsule_viscosity_node_position (IBNode* input_node) {
  IBNode* node = input_node;
  coord position = {0};
  foreach_dimension ()
    position.x = ibval (npos.x);
  return position;
}

static coord ib_capsule_viscosity_triangle_area_normal (IBMesh* mesh,
                                                        int node_ids[3]) {
  coord p[3];
  for (int j = 0; j < 3; j++)
    p[j] = ib_capsule_viscosity_node_position (
      mesh->nodes.ptrs[(size_t) node_ids[j]]);

  coord e10 = {0}, e20 = {0}, area_normal = {0};
  foreach_dimension () {
    e10.x = p[1].x - p[0].x;
    e20.x = p[2].x - p[0].x;
  }

#if dimension == 3
  area_normal.x = 0.5 * (e10.y * e20.z - e10.z * e20.y);
  area_normal.y = 0.5 * (e10.z * e20.x - e10.x * e20.z);
  area_normal.z = 0.5 * (e10.x * e20.y - e10.y * e20.x);
#endif

  return area_normal;
}

void ib_capsule_construct_viscosity_indicator (scalar indicator,
                                               ib_capsule_t capsule,
                                               int mesh_id) {
  assert (capsule);
  assert (mesh_id >= 0 && mesh_id < ibmm.nm);
  IBMesh* mesh = &ibmm.meshes[mesh_id];

  foreach () {
    if (cm[] > 1.e-20) {
      ib_viscosity_divergence[] = 0.;
      foreach_dimension ()
        ib_viscosity_grid_gradient.x[] = 0.;
    }
  }

  const int triangle_count = ib_capsule_get_triangle_count (capsule);
  for (int tid = 0; tid < triangle_count; tid++) {
    int node_ids[3];
    ib_capsule_get_triangle_node_ids (capsule, tid, node_ids);

    if (node_ids[0] < 0 || node_ids[1] < 0 || node_ids[2] < 0)
      continue;
    if ((size_t) node_ids[0] >= mesh->nodes.size ||
        (size_t) node_ids[1] >= mesh->nodes.size ||
        (size_t) node_ids[2] >= mesh->nodes.size)
      continue;

    const coord area_normal =
      ib_capsule_viscosity_triangle_area_normal (mesh, node_ids);

    for (int j = 0; j < 3; j++) {
      IBNode* node = mesh->nodes.ptrs[(size_t) node_ids[j]];
      peskin_cosine_kernel_spread_dimensionless (node) {
        foreach_dimension ()
          ib_viscosity_grid_gradient.x[] -=
            weight * area_normal.x / (3. * dv ());
      }
    }
  }

  foreach_dimension ()
    ib_viscosity_grid_gradient.x.dirty = true;
  boundary ((scalar*) {ib_viscosity_grid_gradient});

  foreach () {
    if (cm[] > 1.e-20) {
      ib_viscosity_divergence[] = 0.;
      foreach_dimension ()
        ib_viscosity_divergence[] +=
          (ib_viscosity_grid_gradient.x[1] -
           ib_viscosity_grid_gradient.x[-1]) /
          (2. * Delta);
    }
  }

  poisson (indicator, ib_viscosity_divergence);

  foreach ()
    if (cm[] > 1.e-20)
      indicator[] = clamp (indicator[], 0., 1.);

  indicator.dirty = true;
  boundary ({indicator});
}

void ib_capsule_set_viscosity_from_indicator (face vector viscosity,
                                              scalar indicator,
                                              double outside_viscosity,
                                              double inside_viscosity) {
  foreach_face ()
    if (fm.x[] > 1.e-20)
      viscosity.x[] =
        (outside_viscosity +
         (inside_viscosity - outside_viscosity) *
           0.5 * (indicator[] + indicator[-1])) *
        fm.x[];
}
