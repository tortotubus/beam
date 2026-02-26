
#include "c/models/ibm/IBMesh.h"
#include <stdlib.h>

extern "C"
{
  /**
   * @memberof ib_mesh_t
   */
  void ib_mesh_free(ib_mesh_t* mesh)
  {
    free(mesh->position);
    free(mesh->velocity);
    free(mesh->forces);
  }
} // extern "C"
