#pragma once

#include "library/ibm/IBMeshManager.h"
#include "library/io/output-dump.h"
#include <errno.h>

struct IBDumpHeader {
  int version, nscalars, nmeshes, nnodes;
};
static const int ib_dump_version = 2;

static IBscalar* ib_dump_list (IBscalar* list);
static bool ib_dump_header (FILE* fp, struct IBDumpHeader* header, IBscalar* list);
static IBscalar ib_scalar_from_name (const char* name);
static int ib_checkpoint_dump (const char* path, void* ctx);
static int ib_checkpoint_restore (const char* path, void* ctx);

void ib_dump (const char* file, IBscalar* list, FILE* fp, bool unbuffered, bool zero);
bool ib_restore (const char* file, IBscalar* list, FILE* fp);

event defaults (i = 0) {
  checkpointer_register ((Checkpointer) {
    .filename = ".ib", .required = true, .dump_phase = CKPT_PHASE_POST_DUMP,
    .dump = ib_checkpoint_dump, .restore_phase = CKPT_PHASE_POST_RESTORE,
    .restore = ib_checkpoint_restore, .ctx = NULL});
}

static int ib_checkpoint_global_ok (int ok) {
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
  return ok;
}

/* Rows are ordered by mesh, then logical node id, never by pool offset.
   Each v2 row stores depth, xyz and the named scalar values. */
static int ib_checkpoint_write (const char* file, IBscalar* list, FILE* stream,
                                 bool unbuffered) {
  IBscalar* fields = ib_dump_list (list);
  fields = iblist_add (fields, nweight);
  int ns = iblist_len (fields), nn = 0;
  foreach_ibmesh () nn += mesh->nnode_global;
  int stride = ns + 4;
  double* data = (double*) calloc ((size_t) nn * stride + 1, sizeof (double));
  int* sources = (int*) malloc (((size_t) nn + 1) * sizeof (int));
  for (int j = 0; j < nn; ++j) sources[j] = INT_MAX;
#if _MPI
  ibmm.dirty = true;
  ibmeshmanager_update_pid ();
#endif
  int offset = 0;
  foreach_ibmesh () {
    for (size_t ni = 0; ni < mesh->nodes.size; ++ni) {
      IBNode* node = mesh->nodes.ptrs[ni];
      /* Prefer center owner; a replica can supply nodes outside the domain. */
      sources[offset + node->node_lid] =
        node->pid == pid () ? pid () : npe () + pid ();
    }
    offset += mesh->nnode_global;
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, sources, nn, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
  int ok = 1;
  for (int j = 0; j < nn; ++j) if (sources[j] == INT_MAX) ok = 0;
  if (ok) {
    offset = 0;
    foreach_ibmesh () {
      for (size_t ni = 0; ni < mesh->nodes.size; ++ni) {
        IBNode* node = mesh->nodes.ptrs[ni];
        int j = offset + node->node_lid;
        if (sources[j] % npe () != pid ()) continue;
        double* row = data + (size_t) j * stride;
        row[0] = node->depth;
        row[1] = ibval (npos.x);
#if dimension > 1
        row[2] = ibval (npos.y);
#endif
#if dimension > 2
        row[3] = ibval (npos.z);
#endif
        int si = 4;
        foreach_ibscalar (fields) row[si++] = ibval (s);
      }
      offset += mesh->nnode_global;
    }
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, data, nn * stride, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (pid () == 0) {
      FILE* fp = stream;
      char* temporary = NULL;
      if (!fp && file) {
        temporary = (char*) malloc (strlen (file) + 2);
        sprintf (temporary, "%s%s", file, unbuffered ? "" : "~");
        fp = fopen (temporary, "wb");
      }
      if (!fp) ok = 0;
      if (ok) {
        struct IBDumpHeader header = {ib_dump_version, ns, ibmm.nm, nn};
        if (!ib_dump_header (fp, &header, fields)) ok = 0;
        foreach_ibmesh () {
          int descriptor[2] = {mesh->gid, mesh->nnode_global};
          if (fwrite (descriptor, sizeof (int), 2, fp) != 2) ok = 0;
        }
        if (fwrite (data, sizeof (double), (size_t) nn * stride, fp) != (size_t) nn * stride)
          ok = 0;
        if (fflush (fp) != 0 || ferror (fp)) ok = 0;
        if (!stream && fclose (fp) != 0) ok = 0;
        if (ok && temporary && !unbuffered && rename (temporary, file) != 0) ok = 0;
      }
      free (temporary);
    }
  }
  free (data); free (sources); free (fields);
  return ib_checkpoint_global_ok (ok) ? 0 : -1;
}

trace void ib_dump (const char* file = "dump.ib", IBscalar* list = iball,
                    FILE* fp = NULL, bool unbuffered = false, bool zero = true) {
  NOT_UNUSED (zero);
  if (ib_checkpoint_write (file, list, fp, unbuffered) != 0) {
    fprintf (stderr, "ib_dump(): unable to write complete logical node state\n");
    exit (1);
  }
}

trace bool ib_restore (const char* file = "dump.ib", IBscalar* list = NULL,
                       FILE* fp = NULL) {
  int nn = 0;
  foreach_ibmesh () nn += mesh->nnode_global;
  struct IBDumpHeader header = {0};
  int ok = 1;
  FILE* input = NULL;
  if (pid () == 0) {
    input = fp ? fp : (file ? fopen (file, "rb") : NULL);
    if (!input) ok = 0;
    if (ok && fread (&header, sizeof (header), 1, input) != 1) ok = 0;
    if (ok && ((header.version != 1 && header.version != 2) ||
               header.nmeshes != ibmm.nm || header.nnodes != nn ||
               header.nscalars < 0 || header.nscalars > 65536)) ok = 0;
  }
#if _MPI
  MPI_Bcast (&header, sizeof (header), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  if (!ib_checkpoint_global_ok (ok)) {
    if (pid () == 0 && input && !fp) fclose (input);
    return false;
  }
  int ns = header.nscalars, stride = ns + 4;
  int* mapping = (int*) malloc (((size_t) ns + 1) * sizeof (int));
  double* data = (double*) calloc ((size_t) nn * stride + 1, sizeof (double));
  IBscalar* allowed = list ? ib_dump_list (list) : NULL;
  if (pid () == 0) {
    for (int si = 0; si < ns && ok; ++si) {
      unsigned len = 0;
      if (fread (&len, sizeof (len), 1, input) != 1 || len > 4096) { ok = 0; break; }
      char* name = (char*) malloc ((size_t) len + 1);
      if (fread (name, 1, len, input) != len) ok = 0;
      name[len] = 0;
      IBscalar field = ib_scalar_from_name (name);
      mapping[si] = field.i;
      free (name);
    }
    if (header.version == 2) {
      foreach_ibmesh () {
        int descriptor[2];
        if (!ok || fread (descriptor, sizeof (int), 2, input) != 2 ||
            descriptor[0] != mesh->gid || descriptor[1] != mesh->nnode_global)
          ok = 0;
      }
      if (ok && fread (data, sizeof (double), (size_t) nn * stride, input) != (size_t) nn * stride)
        ok = 0;
    }
    else {
      /* v1 dense files used registration order and ascending local node ids. */
      int offset = 0;
      foreach_ibmesh () {
        for (int j = 0; j < mesh->nnode_global && ok; ++j) {
          double* row = data + (size_t) (offset + j) * stride;
          row[0] = mesh->depth;
          if (fread (row + 4, sizeof (double), ns, input) != (size_t) ns) ok = 0;
          int found = 0;
          for (int si = 0; si < ns; ++si) {
            if (mapping[si] == npos.x.i) { row[1] = row[4 + si]; found++; }
#if dimension > 1
            if (mapping[si] == npos.y.i) { row[2] = row[4 + si]; found++; }
#endif
#if dimension > 2
            if (mapping[si] == npos.z.i) { row[3] = row[4 + si]; found++; }
#endif
          }
          if (found != dimension) ok = 0;
        }
        offset += mesh->nnode_global;
      }
    }
    if (!fp) fclose (input);
  }
  if (!ib_checkpoint_global_ok (ok)) {
    free (data); free (mapping); free (allowed);
    return false;
  }
#if _MPI
  MPI_Bcast (mapping, ns, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (data, nn * stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  /* Validate before changing any node storage. */
  for (int j = 0; j < nn; ++j) {
    double* row = data + (size_t) j * stride;
    if (!isfinite (row[0]) || row[0] < 0 || row[0] > 30 || row[0] != (int) row[0]) ok = 0;
    for (int k = 1; k < stride; ++k) if (!isfinite (row[k])) ok = 0;
  }
  if (ok) {
    int offset = 0;
    foreach_ibmesh () {
      int epoch = 1;
#if _MPI
      epoch = ++ibmeshmanager_sparse_epoch;
#endif
      for (int j = 0; j < mesh->nnode_global; ++j) {
        double* row = data + (size_t) (offset + j) * stride;
        coord pos = {row[1], row[2], row[3]};
#if _MPI
        if (mesh->sparse_managed &&
            !ibmeshmanager_has_local_kernel_support (pos, (int) row[0]) &&
            !ibmeshmanager_node_center_is_local (pos, (int) row[0])) continue;
#endif
        IBNode* node = ibmesh_get_or_create_node_by_lid (mesh, &ibmm.pool, j);
        node->depth = (int) row[0];
        node->touched_epoch = epoch;
        ibval (npos.x) = pos.x;
#if dimension > 1
        ibval (npos.y) = pos.y;
#endif
#if dimension > 2
        ibval (npos.z) = pos.z;
#endif
        for (int si = 0; si < ns; ++si) {
          IBscalar s = {mapping[si]};
          if (s.i >= 0 && (list == NULL || list == iball || iblist_lookup (allowed, s)))
            ibval (s) = row[4 + si];
        }
      }
      ibmesh_prune_untouched_nodes (mesh, &ibmm.pool, epoch);
      offset += mesh->nnode_global;
    }
    ibmm.dirty = true;
#if _MPI
    ibmeshmanager_update_pid ();
#endif
  }
  free (data); free (mapping); free (allowed);
  return ok != 0;
}

static int ib_checkpoint_dump (const char* path, void* ctx) {
  NOT_UNUSED (ctx);
  return ib_checkpoint_write (path, iball, NULL, false);
}
static int ib_checkpoint_restore (const char* path, void* ctx) {
  NOT_UNUSED (ctx);
  if (ib_restore (path, iball, NULL)) return 0;
  return -1;
}

static IBscalar* ib_dump_list (IBscalar* lista) {
  IBscalar* list = NULL;
  IBscalar* listb = iblist_copy (lista ? lista : iball);

  foreach_ibscalar (listb) {
    if (!ibnodump (s))
      list = iblist_add (list, s);
  }

  free (listb);
  return list;
}

/**
 * @brief
 */
static bool
ib_dump_header (FILE* fp, struct IBDumpHeader* header, IBscalar* list) {
  if (fwrite (header, sizeof (struct IBDumpHeader), 1, fp) < 1) {
    perror ("ib_dump(): error while writing header");
    return false;
  }

  foreach_ibscalar (list) {
    unsigned len = strlen (ibname (s));
    if (fwrite (&len, sizeof (unsigned), 1, fp) < 1) {
      perror ("ib_dump(): error while writing len");
      return false;
    }
    if (fwrite (ibname (s), sizeof (char), len, fp) < len) {
      perror ("ib_dump(): error while writing field name");
      return false;
    }
  }
  return true;
}

/**
 * @brief
 */
static IBscalar ib_scalar_from_name (const char* name) {
  for (size_t i = 0; i < _ibattribute_len; i++) {
    if (_ibattribute[i].name && !strcmp (_ibattribute[i].name, name))
      return (IBscalar) {.i = (int) i};
  }

  return (IBscalar) {.i = -1};
}
