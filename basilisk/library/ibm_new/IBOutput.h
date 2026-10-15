
// #include "output.h"
// #include <limits.h>

#include "io/output-dump.h"

// ============================================================================
// Type declarations
// ============================================================================

struct IBDumpHeader {
  int version;
  int nscalars;
  int nmeshes;
  int nnodes;
};

// ============================================================================
// Globals
// ============================================================================

static const int ib_dump_version = 1;

// ============================================================================
// Function declarations
// ============================================================================

void ib_dump (
  const char* file, IBscalar* list, FILE* fp, bool unbuffered, bool zero);
bool ib_restore (const char* file, IBscalar* list, FILE* fp);

static int ib_checkpoint_dump (const char* path, void* ctx);
static int ib_checkpoint_restore (const char* path, void* ctx);

static bool ib_output_has_sparse_managed_meshes (void);
static IBscalar* ib_dump_list (IBscalar* lista);
static void
ib_dump_header (FILE* fp, struct IBDumpHeader* header, IBscalar* list);

static IBscalar ib_scalar_from_name (const char* name);
static IBscalar*
ib_restore_header_list (FILE* fp, struct IBDumpHeader* header, IBscalar* list);

// ============================================================================
// Events
// ============================================================================

event defaults (i = 0) {
  checkpointer_register (
    (Checkpointer) {.filename = ".ib",
                    .dump_phase = CKPT_PHASE_POST_DUMP,
                    .dump = ib_checkpoint_dump,
                    .restore_phase = CKPT_PHASE_POST_RESTORE,
                    .restore = ib_checkpoint_restore,
                    .ctx = NULL});
}

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 */
static int ib_checkpoint_dump (const char* path, void* ctx) {
  (void) ctx;
  if (ib_output_has_sparse_managed_meshes ())
    return 0;
  ib_dump (path, iball, NULL, false, true);
  return 0;
}

/**
 * @brief
 */
static int ib_checkpoint_restore (const char* path, void* ctx) {
  (void) ctx;
  if (ib_output_has_sparse_managed_meshes ())
    return 0;
  return ib_restore (path, iball, NULL) ? 0 : -1;
}

static bool ib_output_has_sparse_managed_meshes (void) {
  foreach_ibmesh ()
    if (mesh->sparse_managed)
      return true;
  return false;
}

/**
 * @brief
 */
trace void ib_dump (const char* file = "dump.ib",
                    IBscalar* list = iball,
                    FILE* fp = NULL,
                    bool unbuffered = false,
                    bool zero = true) {
  if (ib_output_has_sparse_managed_meshes ())
    return;

  IBscalar* dlist = ib_dump_list (list);
  NOT_UNUSED (zero);
  double* mpi_dump_vals = NULL;
  int mpi_dump_nscalars = iblist_len (dlist);

#if _MPI
  int nscalars = mpi_dump_nscalars;
  if (nscalars > 0 && ibmm.pool.active.size > 0) {
    int nowned = 0;
    foreach_ibnode () {
      if (node->pid == pid ())
        nowned++;
    }

    int* send_ids =
      nowned > 0 ? (int*) malloc ((size_t) nowned * sizeof (int)) : NULL;
    double* send_vals =
      nowned > 0 ? (double*) malloc ((size_t) nowned * (size_t) nscalars *
                                     sizeof (double))
                 : NULL;
    assert (nowned == 0 || (send_ids && send_vals));

    int ni = 0;
    foreach_ibnode () {
      if (node->pid != pid ())
        continue;

      send_ids[ni] = (int) node_id;

      int si = 0;
      foreach_ibscalar (dlist) {
        send_vals[(size_t) ni * (size_t) nscalars + (size_t) si++] = ibval (s);
      }
      ni++;
    }

    int* recv_counts =
      pid () == 0 ? (int*) malloc ((size_t) npe () * sizeof (int)) : NULL;
    MPI_Gather (
      &nowned, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* recv_displs = NULL;
    int* recv_counts_vals = NULL;
    int* recv_displs_vals = NULL;
    int* recv_ids = NULL;
    double* recv_vals = NULL;

    if (pid () == 0) {
      recv_displs = (int*) malloc ((size_t) npe () * sizeof (int));
      recv_counts_vals = (int*) malloc ((size_t) npe () * sizeof (int));
      recv_displs_vals = (int*) malloc ((size_t) npe () * sizeof (int));
      assert (recv_counts && recv_displs && recv_counts_vals &&
              recv_displs_vals);

      int total = 0;
      int total_vals = 0;
      for (int peer = 0; peer < npe (); peer++) {
        recv_displs[peer] = total;
        recv_counts_vals[peer] = recv_counts[peer] * nscalars;
        recv_displs_vals[peer] = total_vals;
        total += recv_counts[peer];
        total_vals += recv_counts_vals[peer];
      }

      recv_ids =
        total > 0 ? (int*) malloc ((size_t) total * sizeof (int)) : NULL;
      recv_vals = total_vals > 0
                    ? (double*) malloc ((size_t) total_vals * sizeof (double))
                    : NULL;
      assert ((total == 0 || recv_ids) && (total_vals == 0 || recv_vals));
    }

    MPI_Gatherv (send_ids,
                 nowned,
                 MPI_INT,
                 recv_ids,
                 recv_counts,
                 recv_displs,
                 MPI_INT,
                 0,
                 MPI_COMM_WORLD);

    MPI_Gatherv (send_vals,
                 nowned * nscalars,
                 MPI_DOUBLE,
                 recv_vals,
                 recv_counts_vals,
                 recv_displs_vals,
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

    if (pid () == 0) {
      int total = 0;
      for (int peer = 0; peer < npe (); peer++)
        total += recv_counts[peer];

      mpi_dump_vals = (double*) malloc ((size_t) ibmm.pool.active.size *
                                        (size_t) nscalars * sizeof (double));
      assert (mpi_dump_vals);

      foreach_ibnode () {
        int si = 0;
        foreach_ibscalar (dlist) {
          mpi_dump_vals[(size_t) node_id * (size_t) nscalars + (size_t) si++] =
            ibval (s);
        }
      }

      for (int i = 0; i < total; i++) {
        int node_id_remote = recv_ids[i];
        assert (node_id_remote >= 0 &&
                node_id_remote < (int) ibmm.pool.active.size);
        for (int si = 0; si < nscalars; si++)
          mpi_dump_vals[(size_t) node_id_remote * (size_t) nscalars +
                        (size_t) si] =
            recv_vals[(size_t) i * (size_t) nscalars + (size_t) si];
      }
    }

    free (send_ids);
    free (send_vals);
    free (recv_counts);
    free (recv_displs);
    free (recv_counts_vals);
    free (recv_displs_vals);
    free (recv_ids);
    free (recv_vals);
  }
#endif

  if (pid () == 0) {
    char* name = NULL;
    if (!fp) {
      name = (char*) malloc (strlen (file) + 2);
      strcpy (name, file);
      if (!unbuffered)
        strcat (name, "~");
      if ((fp = fopen (name, "wb")) == NULL) {
        perror (name);
        exit (1);
      }
    }
    assert (fp);
    struct IBDumpHeader header = {.version = ib_dump_version,
                                  .nscalars = iblist_len (dlist),
                                  .nmeshes = ibmm.nm,
                                  .nnodes = (int) ibmm.pool.active.size};

    ib_dump_header (fp, &header, dlist);

    foreach_ibnode () {
      int si = 0;
      foreach_ibscalar (dlist) {
        double val = mpi_dump_vals
                       ? mpi_dump_vals[(size_t) node_id *
                                         (size_t) mpi_dump_nscalars +
                                       (size_t) si]
                       : ibval (s);
        si++;
        if (fwrite (&val, sizeof (double), 1, fp) < 1) {
          perror ("ib_dump(): error while writing scalars");
          exit (1);
        }
      }
    }

    if (file) {
      fclose (fp);
      if (!unbuffered)
        rename (name, file);
      free (name);
    }
  }

  free (mpi_dump_vals);
  free (dlist);
}

/**
 * @brief
 */
trace bool ib_restore (const char* file = "dump.ib",
                       IBscalar* list = NULL,
                       FILE* fp = NULL) {

  if (ib_output_has_sparse_managed_meshes ())
    return true;

  if (!fp && (fp = fopen (file, "rb")) == NULL)
    return false;
  assert (fp);

  struct IBDumpHeader header = {0};
  if (fread (&header, sizeof (header), 1, fp) < 1) {
    fprintf (ferr, "ib_restore(): error: expecting header\n");
    exit (1);
  }

  if (header.version != ib_dump_version) {
    fprintf (ferr,
             "ib_restore(): error: file version mismatch: %d (file) != %d "
             "(code)\n",
             header.version,
             ib_dump_version);
    exit (1);
  }

  if (header.nmeshes != ibmm.nm) {
    fprintf (ferr,
             "ib_restore(): error: mesh count mismatch: %d (file) != %d "
             "(code)\n",
             header.nmeshes,
             ibmm.nm);
    exit (1);
  }

  if (header.nnodes != (int) ibmm.pool.active.size) {
    fprintf (ferr,
             "ib_restore(): error: node count mismatch: %d (file) != %d "
             "(code)\n",
             header.nnodes,
             (int) ibmm.pool.active.size);
    exit (1);
  }

  IBscalar* slist = ib_restore_header_list (fp, &header, list);

  foreach_ibnode () {
    foreach_ibscalar (slist) {
      double val = 0.;
      if (fread (&val, sizeof (double), 1, fp) < 1) {
        fprintf (ferr, "ib_restore(): error: expecting scalar\n");
        exit (1);
      }
      if (s.i != INT_MAX)
        ibval (s) = val;
    }
  }
  if (file)
    fclose (fp);

  free (slist);

#if _MPI
  ibmm.dirty = true;
  ibmeshmanager_update_pid();
#endif 

  return true;
}

/**
 * @brief
 */
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
static void
ib_dump_header (FILE* fp, struct IBDumpHeader* header, IBscalar* list) {
  if (fwrite (header, sizeof (struct IBDumpHeader), 1, fp) < 1) {
    perror ("ib_dump(): error while writing header");
    exit (1);
  }

  foreach_ibscalar (list) {
    unsigned len = strlen (ibname (s));
    if (fwrite (&len, sizeof (unsigned), 1, fp) < 1) {
      perror ("ib_dump(): error while writing len");
      exit (1);
    }
    if (fwrite (ibname (s), sizeof (char), len, fp) < len) {
      perror ("ib_dump(): error while writing field name");
      exit (1);
    }
  }
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

/**
 * @brief
 */
static IBscalar*
ib_restore_header_list (FILE* fp, struct IBDumpHeader* header, IBscalar* list) {
  IBscalar* input = NULL;
  IBscalar* allowed = list ? ib_dump_list (list) : NULL;
  bool restore_all = (list == NULL || list == iball);

  for (int i = 0; i < header->nscalars; i++) {
    unsigned len = 0;
    if (fread (&len, sizeof (unsigned), 1, fp) < 1) {
      fprintf (ferr, "ib_restore(): error: expecting len\n");
      exit (1);
    }

    char name[len + 1];
    if (fread (name, sizeof (char), len, fp) < len) {
      fprintf (ferr, "ib_restore(): error: expecting field name\n");
      exit (1);
    }
    name[len] = '\0';

    IBscalar s = ib_scalar_from_name (name);
    if (s.i >= 0 && (restore_all || iblist_lookup (allowed, s)))
      input = iblist_append (input, s);
    else
      input = iblist_append (input, (IBscalar) {.i = INT_MAX});
  }

  free (allowed);
  return input;
}
