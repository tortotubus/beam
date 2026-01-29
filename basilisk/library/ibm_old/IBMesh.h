#include "interface/ibm/immersedboundary.h"

#include "library/ibm/IBNode.h"
#include "library/ibm/IBNodeList.h"

#include <stdio.h>

/* ====================================================================================================================
 * IBMesh
 * ====================================================================================================================
 */

#define RHO_F 1

typedef struct {
  ib_force_structure_model_t model;
  // int* model;
  int* node_ids;
  int nn;
  int pid;
} IBMesh;

/*
 * Macro defintions
 */

macro foreach_ibnode_mesh (IBMesh* mesh) {
  for (int _nii = 0; _nii < mesh->nn; _nii++) {
    int node_index = mesh->node_ids[_nii];
    IBNode* node = ibnode_list_get_node (node_index);
    NOT_UNUSED (node);
    // clang-format off
    {...}
    // clang-format on
  }
}

/*
 * Functions
 */

/*!
 * @brief Default constructor for an IBMesh
 *
 * @param mesh Mesh to initialize
 *
 * @memberof IBMesh
 */
void ib_mesh_init (IBMesh* mesh) {
  mesh->model = NULL;
  mesh->node_ids = NULL;
  mesh->nn = 0;
  mesh->pid = -1;
}

/*!
 * @brief Destructor for the IBMesh struct
 * @memberof IBMesh
 */
void ib_mesh_free (IBMesh* mesh) {
  if (mesh->node_ids) {
    free (mesh->node_ids);
    mesh->node_ids = NULL;
    mesh->model = NULL;
    mesh->nn = 0;
    mesh->pid = -1;
  }
}

/*!
 * @brief Construct the IBMesh with a number of IBNodes
 *
 * @param mesh The IBMesh to initialize
 * @param nn The number of nodes to initialize
 *
 * @memberof IBMesh
 */
void ib_mesh_init_nn (IBMesh* mesh, int nn) {
  assert (mesh->node_ids == NULL);

  mesh->node_ids = (int*) calloc (nn, sizeof (int));
  mesh->model = NULL;
  mesh->nn = nn;
  mesh->pid = -1;

  for (int nii = 0; nii < nn; nii++) {
    IBNode new_node;
    ibnode_init (&new_node);
    int node_index = ibnode_list_add_node (new_node);
    mesh->node_ids[nii] = node_index;
  }
}

/*!
 * @brief Set the model for an IBMesh
 *
 * @param mesh The IBMesh
 * @param model The pointer to the force-coupled structural model
 *
 * @memberof IBMesh
 */
// void ib_mesh_model (IBMesh* mesh, ib_force_structure_model_t model) {
// #if DEBUG
//   printf ("ib_mesh_model\n");
// #endif
//   // Free the previous mesh, if any
//   ib_mesh_free (mesh);

//   // Set the model pointer
//   mesh->model = model;

//   // Get the number of nodes from the C++ model
//   int nn = ib_force_structure_model_get_number_of_nodes (model);

//   // Initialize the basilisk C mesh
//   ib_mesh_init_nn (mesh, nn);

//   // Get the current-time mesh from the C++ model
//   ib_structure_mesh_t cpp_mesh = ib_force_structure_model_get_current (model);

//   foreach_ibnode (mesh) {
//     foreach_dimension () {
//       node->lagpos.x = cpp_mesh.position[node_index].x;
//       node->lagvel.x = cpp_mesh.velocity[node_index].x;
//     }
//   }

//   // Free the C++ mesh
//   ib_structure_mesh_free (&cpp_mesh);
// }

/*!
 * @brief Helper function to compute the centroid of a given IBMesh
 *
 * @param mesh IBMesh to compute the centroid of
 *
 * @memberof IBMesh
 */
coord ib_mesh_get_centroid (IBMesh* mesh) {
  coord centroid = {.x = 0, .y = 0, .z = 0};

  foreach_ibnode_mesh (mesh) {
    foreach_dimension () {
      centroid.x += node->lagpos.x;
    }
  }

  foreach_dimension () {
    if (centroid.x != 0)
      centroid.x /= (double) mesh->nn;
  }

  return centroid;
}

/*!
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 * @param mesh
 *
 * @memberof IBMesh
 */

trace void ib_mesh_advance_lagrangian_mesh (IBMesh* mesh) {
  // vertex_t* force = calloc (mesh->nn, sizeof (vertex_t));

  // foreach_ibnode (mesh) {
  //   foreach_dimension () {
  //     force[_nii].x = -node->force.x * node->weight + node->gravity.x;
  //   }
  // }

  // ib_structure_mesh_t cpp_mesh =
  //   ib_force_structure_model_get_next (mesh->model, force, mesh->nn, dt);

  // free (force);

  // foreach_ibnode (mesh) {
  //   foreach_dimension () {
  //     node->lagpos.x = cpp_mesh.position[node_index].x;
  //     node->lagvel.x = cpp_mesh.velocity[node_index].x;
  //   }
  // }

  // ib_structure_mesh_free (&cpp_mesh);
}

