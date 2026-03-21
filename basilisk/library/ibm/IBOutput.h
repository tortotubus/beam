
#include "output.h"
#include <limits.h>

#include "library/elff/runtime.h"
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

void ib_dump (const char* file, IBscalar* list, FILE* fp, bool unbuffered, bool zero);
bool ib_restore (const char* file, IBscalar* list, FILE* fp);

void ib_and_basilisk_dump(const char * file, scalar * list, FILE* fp, bool unbuffered, bool zero);
bool ib_and_basilisk_restore(const char * file, scalar * list, FILE *fp);

static IBscalar* ib_dump_list (IBscalar* lista);
static void ib_dump_header (FILE* fp, struct IBDumpHeader* header, IBscalar* list);

static IBscalar ib_scalar_from_name (const char* name);
static IBscalar* ib_restore_header_list (FILE* fp, struct IBDumpHeader* header, IBscalar* list);

// ============================================================================
// Macros
// ============================================================================

@define dump_function(fname,slist) ib_and_basilisk_dump((fname),(slist),NULL,false,true)
@define restore_function(fname,slist) ib_and_basilisk_restore((fname),(slist),NULL)

#include "library/io/output-dump.h"

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 */
void ib_and_basilisk_dump(const char * file = "dump", scalar * list = all, FILE * fp = NULL, bool unbuffered = false, bool zero = true) {
  char * ibname = (char *) malloc(strlen(file) + 4);
  strcpy(ibname,file);
  strcat(ibname,".ib");

  char * elffname = (char *) malloc(strlen(file) + 6);
  strcpy(elffname,file);
  strcat(elffname,".elff");

  if (elff_dump(elffname) != 0) {
    fprintf (ferr, "ib_and_basilisk_dump(): error: failed to write ELFF checkpoint '%s'\n",
             elffname);
    exit (1);
  }
  ib_dump(ibname, iball, NULL, unbuffered, zero);
  dump(file,list,fp,unbuffered,zero);

  free(elffname);
  free(ibname);
}

/**
 * @brief
 */
bool ib_and_basilisk_restore(const char * file = "dump", scalar * list = all, FILE * fp = NULL) {
  char * ibname = (char *) malloc(strlen(file) + 4);
  strcpy(ibname,file);
  strcat(ibname,".ib");

  char * elffname = (char *) malloc(strlen(file) + 6);
  strcpy(elffname,file);
  strcat(elffname,".elff");

  if (elff_restore(elffname) != 0) {
    fprintf (ferr, "ib_and_basilisk_restore(): error: failed to restore ELFF checkpoint '%s'\n",
             elffname);
    exit (1);
  }
  if (!ib_restore(ibname, iball, NULL)) {
    fprintf (ferr, "ib_and_basilisk_restore(): error: failed to restore IB checkpoint '%s'\n",
             ibname);
    exit (1);
  }
  bool ok = restore(file,list,fp);

  free(elffname);
  free(ibname);
  return ok;
}


#if !_MPI
/**
 * @brief
 */
trace void ib_dump (const char* file = "ibdump",
                    IBscalar* list = iball,
                    FILE* fp = NULL,
                    bool unbuffered = false,
                    bool zero = true) {
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

  IBscalar* dlist = ib_dump_list (list);
  struct IBDumpHeader header = {.version = ib_dump_version,
                                .nscalars = iblist_len (dlist),
                                .nmeshes = ibmm.nm,
                                .nnodes = (int) ibmm.pool.active.size};

  ib_dump_header (fp, &header, dlist);

  foreach_ibnode () {
    foreach_ibscalar (dlist) {
      double val = ibval (s);
      if (fwrite (&val, sizeof (double), 1, fp) < 1) {
        perror ("ib_dump(): error while writing scalars");
        exit (1);
      }
    }
  }

  free (dlist);

  if (file) {
    fclose (fp);
    if (!unbuffered)
      rename (name, file);
    free (name);
  }
}

/**
 * @brief
 */
trace bool ib_restore (const char* file = "ibdump",
                       IBscalar* list = NULL,
                       FILE* fp = NULL) {
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

  free (slist);

  if (file)
    fclose (fp);

  return true;
}
#else // _MPI
#endif


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
