#include <stdio.h>

#include "library/ibm/IBNode.h"
#include "library/ibm/IBMesh.h"
#include "library/ibm/IBMeshManager.h"

void output_ibnodes (const char* basename, int iter, double time) {

  int total_points = 0;
  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      total_points++;
    }
  }

  // if (total_points == 0) {
  //   return;
  // }

  // Name the file
  char fname[128];
  sprintf (fname, "%s_pid_%d_%d.vtk", basename, pid (), iter);

  // Open the file
  FILE* fp = fopen (fname, "w");

  if (!fp) {
    char error_msg[256];
    sprintf (error_msg, "error: fopen for %s failed\n", fname);
    fprintf (stderr, error_msg);
    return;
  }

  // Legacy VTK header
  fprintf (fp, "# vtk DataFile Version 3.0\n");
  fprintf (fp, "IB points at step %d\n", iter);
  fprintf (fp, "ASCII\n");
  fprintf (fp, "DATASET POLYDATA\n");
  fprintf (fp, "POINTS %d double\n", total_points);

  // Write point coordinate(s)
  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (fp,
               "%.16g %.16g %.16g\n",
               node->lagpos.x,
               node->lagpos.y,
               node->lagpos.z);
    }
  }

  // Write vertices for each point
  fprintf (fp, "VERTICES %d %d\n", total_points, 2 * total_points);
  for (int i = 0; i < total_points; i++) {
    fprintf (fp, "1 %d\n", i);
  }

  // ===========================
  // POINT DATA (one value per node)
  // ===========================
  fprintf (fp, "\nPOINT_DATA %d\n", total_points);

  // Example 1: scalar “mesh_id” telling which mesh the node belongs to
  // fprintf (fp, "SCALARS mesh_id int 1\n");
  // fprintf (fp, "LOOKUP_TABLE default\n");

  // foreach_ibnode(mesh)() {

  //     // assume mesh has an integer id; if not, use a counter
  //     fprintf (fp, "%d\n", mesh_index);

  // }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS force double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf ( fp, "%.16g %.16g %.16g\n", node->force.x, node->force.y, node->force.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS lagvel double\n");
  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (fp,
               "%.16g %.16g %.16g\n",
               node->lagvel.x,
               node->lagvel.y,
               node->lagvel.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS eulvel double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (fp,
               "%.16g %.16g %.16g\n",
               node->eulvel.x,
               node->eulvel.y,
               node->eulvel.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS rhs double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (
        fp, "%.16g %.16g %.16g\n", node->rhs.x, node->rhs.y, node->rhs.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS res double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (
        fp, "%.16g %.16g %.16g\n", node->res.x, node->res.y, node->res.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS w double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (fp, "%.16g %.16g %.16g\n", node->w.x, node->w.y, node->w.z);
    }
  }

  // Example 2: vector “lag_velocity” at each node
  fprintf (fp, "VECTORS Ay double\n");

  foreach_ibmesh () {
    foreach_ibnode (mesh)  {
      fprintf (fp, "%.16g %.16g %.16g\n", node->Ay.x, node->Ay.y, node->Ay.z);
    }
  }

  fclose (fp);
}