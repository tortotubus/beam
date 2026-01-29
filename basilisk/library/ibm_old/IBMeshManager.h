#include "library/ibm/IBMesh.h"

/*
 * Forward declaration of this header's functions
 */
void ib_mesh_manager_init (int nm);
void ib_mesh_manager_free ();
// void ib_mesh_manager_generate_stencil_caches ();
// void ib_mesh_manager_tag_stencils ();

/*
 * Type definitions
 */

/*!
 * @brief Management struct to hold all IBMeshes in the simulation
 */
typedef struct {
  IBMesh* meshes; /*!< The actual immersed boundary meshes */
  int nm;         /*!< Total nuber of meshes*/
} IBMeshManager;

/*
 * Globals
 */

IBMeshManager ib_mesh_manager;

/*
 * Macro definitions
 */

#define ibmesh(i) &ib_mesh_manager.meshes[i]

macro foreach_ibmesh () {
  for (int mesh_index = 0; mesh_index < ib_mesh_manager.nm; mesh_index++) {
    IBMesh* mesh = &ib_mesh_manager.meshes[mesh_index];
    NOT_UNUSED (mesh);
    // clang-format off
    {...}
    // clang-format on
  }
}

/*
 * Function definitions
 */

/*!
 * @brief Initialize the immersed boundary mesh manager
 *
 * Although the meshes are default constructed and can be changed after calling
 * this, the number of meshes must remain fixed throughout the simulation. Since
 * meshes can be set to active or inactive, it is reccomended to choose the
 * number of meshes as the maximum expected to be needed throughout the life of
 * the simulation.
 *
 * @param nm The number of meshes you plan to have
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_init (int nm) {
#if DEBUG
  printf ("ib_mesh_manager_init\n");
#endif
  if (ib_mesh_manager.meshes != NULL) {
    ib_mesh_manager_free ();
  }

  ib_mesh_manager.meshes = NULL;
  ib_mesh_manager.nm = nm;
  ib_mesh_manager.meshes = (IBMesh*) calloc (nm, sizeof (IBMesh));

  // Default constructor for each of our meshes
  // foreach_ibmesh () {
  //   ib_mesh_init (mesh);
  // }
}

/*!
 * @brief Free all members in the immersed boundary mesh manager, including the
 * array of meshes, their nodes, and stencil nodes.
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_free () {
#if DEBUG
  printf ("ib_mesh_manager_free\n");
#endif

  // Call the IBMesh destructor on each mesh
  foreach_ibmesh () {
    ib_mesh_free (mesh);
  }

  // Free the owned arrays
  free (ib_mesh_manager.meshes);

  ib_mesh_manager.meshes = NULL;
  ib_mesh_manager.nm = 0;
}

/*!
 * @brief For force-coupled meshes, we advance the mesh in time according to
 * it's own kinematics using the previous force at time \f(t^{n}\f)
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_advance_lagrangian_mesh () {
  // foreach_ibmesh () {
  //   if (pid () == mesh->pid)
  //     ib_mesh_advance_lagrangian_mesh (mesh);
  // }
}
