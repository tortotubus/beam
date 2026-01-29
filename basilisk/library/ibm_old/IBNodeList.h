#include "library/ibm/IBNode.h"

#include <stdio.h>

/*
 * IBNodeList
 */

typedef struct {
  IBNode* nodes;
  int size;
  int capacity;
} IBNodeList;

/*
 * IBNodeList globals
 */
IBNodeList ibnode_list;

/*
 * IBNodeList macros
 */
macro foreach_ibnode_all (bool local_only = true) {
  for (int _node_index = 0; _node_index < ibnode_list.size; _node_index++) {
    IBNode* node = &ibnode_list.nodes[_node_index];
    if (local_only) {
      if (node->pid == pid ()) {
        // clang-format off
        {...}
        // clang-format on
      }
    } else {
      // clang-format off
      {...}
      // clang-format on      
    }
  }
}

/*
 * IBNodeList functions
 */

/*!
 *
 */
void ibnode_list_free () {

  for (int ni = 0; ni < ibnode_list.capacity; ni++) {
    ibnode_free (&ibnode_list.nodes[ni]);
  }

  free (ibnode_list.nodes);

  ibnode_list.nodes = NULL;
  ibnode_list.capacity = 0;
  ibnode_list.size = 0;
}

/*!
 *
 */
void ibnode_list_init (int capacity) {
  if (capacity <= 0)
    return;

  // ibnode_list_free();

  ibnode_list.nodes = (IBNode*) malloc (capacity * sizeof (IBNode));

  if (!ibnode_list.nodes) {
    fprintf (stderr, "error: malloc failure in ibnode_list_extend");
  }

  ibnode_list.capacity = capacity;
  ibnode_list.size = 0;

  for (int ni = 0; ni < capacity; ni++) {
    ibnode_init (&ibnode_list.nodes[ni]);
  }
}

/*!
 *
 */
void ibnode_list_extend (int capacity_new) {
  IBNode* nodes_old = ibnode_list.nodes;
  int capacity_old = ibnode_list.capacity;
  int size_old = ibnode_list.size;

  assert (capacity_new >= capacity_old);

  IBNode* nodes_new = (IBNode*) malloc (capacity_new * sizeof (IBNode));

  if (!nodes_new) {
    fprintf (stderr, "error: malloc failure in ibnode_list_extend");
  }

  for (int ni = 0; ni < size_old; ni++) {
    ibnode_copy (&nodes_new[ni], &nodes_old[ni]);
  }

  for (int ni = size_old; ni < capacity_new; ni++) {
    ibnode_init (&nodes_new[ni]);
  }

  ibnode_list.nodes = nodes_new;
  ibnode_list.capacity = capacity_new;
  ibnode_list.size = size_old;

  free (nodes_old);
}

/*!
 *
 * @returns Index where the node is placed
 */
int ibnode_list_add_node (IBNode node) {
  if (ibnode_list.size == ibnode_list.capacity) {
    ibnode_list_extend (ibnode_list.capacity * 2);
  }

  int node_id = ibnode_list.size;
  ibnode_list.nodes[node_id] = node;
  ibnode_list.size++;

  return node_id;
}

/*!
 *
 */
IBNode * ibnode_list_get_node(int index) {
  if (index >= ibnode_list.size) {
    return NULL;
  } else {
    return &ibnode_list.nodes[index];
  }
}

/*!
 * 
 */
void ibnode_list_spread_eulerian_forcing(vector forcing) {
  foreach_ibnode_all() {
    foreach_neighbor_of_coord(node->lagpos, PESKIN_SUPPORT_RADIUS) {
      peskin_kernel(node->lagpos) {
        foreach_dimension() {
          forcing.x[] += peskin_weight * node->force.x * node->weight;
        }
      }
    }
  }
}


/*!
 * 
 */
void ibnode_list_interpolate_eulerian_velocities() {
  foreach_ibnode_all() {
    foreach_neighbor_of_coord(node->lagpos, PESKIN_SUPPORT_RADIUS) {
      peskin_kernel(node->lagpos) {
        foreach_dimension() {
          node->eulvel.x += peskin_weight * u.x[];
        }
      }
    }
  }
}

/*!
 *
 */
void ibnode_list_set_pid() {
  foreach_ibnode_all(local_only = false) {
    Point p = locate (node->lagpos.x, node->lagpos.y, node->lagpos.z);
    if (p.level > -1) {
      node->pid = pid ();
    } else {
      node->pid = -1;
    }
  }
}

/*!
 *
 */
void ibnode_list_generate_stencil_cache() {
  foreach_ibnode_all(local_only = true) {
    node->stencil.n = 0;
    foreach_neighbor_of_coord(lagpos, PESKIN_SUPPORT_RADIUS) {
      if (point.level >= 0 && point.level != grid->maxdepth) {
        fprintf (
          stderr,
          "warning: lagrangian stencil not fully resolved at level %d of "
          "%d.\n",
          point.level,
          grid->maxdepth);
      }

      cache_append (&node->stencil, point, 0);
    }
  }
}

/*!
 * 
 */
void ibnode_list_compute_constraint_rhs() {
  foreach_ibnode_all() {
    foreach_dimension() {
      node->rhs.x = node->lagvel.x - node->eulvel.x;
    }
  }
}

vector ibnode_cg_velocity[];
vector ibnode_cg_lambda[];

/*!
 *
 */
void ibnode_list_matvec_Aw(double dt) {
  foreach() {
    foreach_dimension() {
      ibnode_cg_lambda.x[] = 0;
      ibnode_cg_velocity.x[] = 0;
    }
  }

  // 1. g = J^T w on grid
  foreach_ibnode_all() {
    foreach_cache(node->stncil) {
      peskin_kernel(node->lagpos) {
        foreach_dimension() {
          ibnode_cg_lambda.x[] += peskin_weight * node->w.x * node->weight;
        }
      }
    }
  }

  // #if _MPI
  //   synchronize ((scalar*) {vector ibnode_cg_velocity});
  // #endif

  //2. c = M_f^{-1} g ~ g / (rho_f)
  foreach () {
    foreach_dimension () {
      ibnode_cg_velocity.x[] = ibnode_cg_lambda.x[] / 1.;
    }
  }

  // #if _MPI
  //   synchronize ((scalar*) {vel});
  // #endif

  // 3. d = J c, store in node->Ay temporarily, then scale by dt
  foreach_ibnode_all() {
    foreach_dimension() {
      node->Ay.x = 0.;
    } 
    foreach_cache(node->stencil) {
      peskin_kernel(node->lagpos) {
        foreach_dimension() {
          node->Ay.x += peskin_weight * ibnode_cg_velocity.x[];
        }
      }
    } 
    foreach_dimension () {
      node->Ay.x *= dt;
    }
  }
}

/*!
 *
 */
void ibnode_list_solve_lambda_CG (double dt) {
  const int maxiter = 5000;
  const double tol = 1e-6;

  foreach_ibnode_all() {
    foreach_dimension() {
      node->force.x = 0;
      node->res.x = node->rhs.x;
      node->w.x = node->res.x;
    }
  }

  double rr_old = 0.;
  foreach_ibnode_all() {
    foreach_dimension() {
      rr_old += sq(node->res.x) * node->weight;
    }
  }

#if _MPI
  double rr_old_local = rr_old;
  MPI_Allreduce (
    &rr_old_local, &rr_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (sqrt (rr_old) < tol) {
    return;
  }

  for (int cg_iter = 0; cg_iter < maxiter; cg_iter++) {
    ibnode_list_matvec_Aw(dt);

    double wy = 0.;
    foreach_ibnode_all() {
      foreach_dimension() {
        wy += node->w.x * node->Ay.x * node->weight;
      }
    }

#if _MPI
    double wy_local = wy;
    MPI_Allreduce (&wy_local, &wy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (fabs(wy) < 1e-30) {
      break;
    }

    double alpha = rr_old / wy;
    double rr_new = 0.;

    foreach_ibnode_all() {
      foreach_dimension() {
        node->force.x += alpha * node->w.x;
        node->res.x -= alpha * node->Ay.x;
        rr_new += sq(node->res.x) * node->weight;
      }
    }

#if _MPI
    double rr_new_local = rr_new;
    MPI_Allreduce (
      &rr_new_local, &rr_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    // Check tolerance
    if (sqrt (rr_new) < tol)
      break;

    double beta = rr_new / rr_old;
    rr_old = rr_new;

    // w_{k+1} = r_{k+1} + \beta w_{k}
    foreach_ibnode_all () {
      foreach_dimension () {
        node->w.x = node->res.x + beta * node->w.x;
      }
    }
  }
}

scalar stencils[];

void ibnode_list_tag_stencil() {
  foreach() {
    stencils[] = 0;
  }

  foreach_ibnode_all() {
    foreach_cache(node->stencil) {
      peskin_kernel(node->lagpos) {
        double length = 0;
        foreach_dimension () {
          length += dist.x;
        }
        length = sq (length);
#if dimension == 1
        stencils[] = length / (2. * Delta) * (2. + noise ());
#elif dimension == 2
        stencils[] = length / sq (2. * Delta) * (2. + noise ());
#else
        stencils[] = length / cube (2. * Delta) * (2. + noise ());
#endif
      }
    }
  }
}