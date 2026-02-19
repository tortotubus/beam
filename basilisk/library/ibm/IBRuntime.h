// #include "common.h"

#include "library/ibm/IBNode.h"
#include "library/ibm/IBNodeList.h"
#include "library/ibm/IBMempool.h"
#include "library/ibm/IBMesh.h"
#include "library/ibm/IBMeshManager.h"

// Create a global instance of IBMeshManager
IBMeshManager ib_mesh_manager = {0};

// Create a free function for solver teardown
void ibm_solver_free() { ibmeshmanager_free(&ib_mesh_manager); }

// Register the free function 
free_solver_func_add(ibm_solver_free);