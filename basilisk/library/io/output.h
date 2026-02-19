/**
## *dump()*: Basilisk snapshots

This function (together with *restore()*) can be used to dump/restore
entire simulations.

The arguments and their default values are:

*file*
: the name of the file to write to (mutually exclusive with *fp*). The
default is "dump".

*list*
: a list of scalar fields to write. Default is *all*. 

*fp*
: a file pointer. Default is stdout.

*unbuffered*
: whether to use a file buffer. Default is false.

*zero*
: whether to dump fields which are zero. Default is true.
*/

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =
  // 161020
  170901;

static scalar * dump_list (scalar * lista, bool zero)
{
  scalar * list = is_constant(cm) ? NULL : list_concat ({cm}, NULL);
  // fixme: on GPUs statsf() can change the `all` list, because it
  // allocates new fields to store reductions, which causes a nasty
  // crash...
#if 1
  scalar * listb = list_copy (lista);
#endif
  for (scalar s in listb)
    if (!s.face && !s.nodump && s.i != cm.i) {
      if (zero)
	list = list_add (list, s);
      else {	
	stats ss = statsf (s);
	if (ss.min != 0. || ss.max != 0.)
	  list = list_add (list, s);
      }
    }
  free (listb);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  for (scalar s in list) {
    unsigned len = strlen(s.name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (s.name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !_MPI
trace
void dump (const char * file = "dump",
	   scalar * list = all,
	   FILE * fp = NULL,
	   bool unbuffered = false,
	   bool zero = true)
{
  char * name = NULL;
  if (!fp) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if (!unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  assert (fp);
  
  scalar * dlist = dump_list (list, zero);
  scalar size[];
  scalar * slist = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
			       dump_version };
  int npe = 1;
  foreach_dimension() {
    header.n.x = Dimensions.x;
    npe *= header.n.x;
  }
  header.npe = npe;
  dump_header (fp, &header, slist);
  
  subtree_size (size, false);
#if _GPU
  for (scalar s in slist)
    s.input = 1;
  gpu_cpu_sync (slist, GL_MAP_READ_BIT, __FILE__, LINENO);
#endif // _GPU
  foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    for (scalar s in slist) {
      double val = s[];
      if (fwrite (&val, sizeof(double), 1, fp) < 1) {
	perror ("dump(): error while writing scalars");
	exit (1);
      }
    }
    if (is_leaf(cell))
      continue;
  }
  
  free (slist);
  if (file) {
    fclose (fp);
    if (!unbuffered)
      rename (name, file);
    free (name);
  }
}
#else // _MPI
trace
void dump (const char * file = "dump",
	   scalar * list = all,
	   FILE * fp = NULL,
	   bool unbuffered = false,
	   bool zero = true)
{
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);    
  }

  scalar * dlist = dump_list (list, zero);
  scalar size[];
  scalar * slist = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
			       dump_version };

#if MULTIGRID_MPI
  foreach_dimension()
    header.n.x = Dimensions.x;
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  if (pid() == 0)
    dump_header (fh, &header, slist);
  
  scalar index = {-1};
  
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  for (scalar s in slist)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(s.name);
  long pos = pid() ? 0 : sizeofheader;
  
  subtree_size (size, false);
  
  foreach_cell() {
    // fixme: this won't work when combining MPI and mask()
    if (is_local(cell)) {
      long offset = sizeofheader + index[]*cell_size;
      if (pos != offset) {
	fseek (fh, offset, SEEK_SET);
	pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      for (scalar s in slist) {
	double val = s[];
	fwrite (&val, 1, sizeof(double), fh);
      }
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }

  delete ({index});
  
  free (slist);
  fclose (fh);
  if (!unbuffered && pid() == 0)
    rename (name, file);
}
#endif // _MPI

trace
bool restore (const char * file = "dump",
	      scalar * list = NULL,
	      FILE * fp = NULL)
{
  if (!fp && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  struct DumpHeader header = {0};
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }

#if TREE
  init_grid (1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
#else // multigrid
#if MULTIGRID_MPI
  if (header.npe != npe()) {
    fprintf (ferr,
	     "restore(): error: the number of processes don't match:"
	     " %d != %d\n",
	     header.npe, npe());
    exit (1);
  }
#endif // MULTIGRID_MPI
  dimensions (header.n.x, header.n.y, header.n.z);
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);
#endif // multigrid

  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all, true);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
	       "restore(): error: the list lengths don't match: "
	       "%ld (file) != %d (code)\n",
	       header.len - 1, list_len (slist));
      exit (1);
    }
  }
  else { // header.version != 161020
    if (header.version != dump_version) {
      fprintf (ferr,
	       "restore(): error: file version mismatch: "
	       "%d (file) != %d (code)\n",
	       header.version, dump_version);
      exit (1);
    }
    
    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting len\n");
	exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting s.name\n");
	exit (1);
      }
      name[len] = '\0';

      if (i > 0) { // skip subtree size
	bool found = false;
	for (scalar s in slist)
	  if (!strcmp (s.name, name)) {
	    input = list_append (input, s);
	    found = true; break;
	  }
	if (!found) {
	  if (restore_all) {
	    scalar s = new scalar;
	    free (s.name);
	    s.name = strdup (name);
	    input = list_append (input, s);
	  }
	  else
	    input = list_append (input, (scalar){INT_MAX});
	}
      }
    }
    free (slist);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }

#if MULTIGRID_MPI
  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << dimension*(header.depth + 1)) - 1)/
    ((1 << dimension) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }
#endif // MULTIGRID_MPI
  
  scalar * listm = is_constant(cm) ? NULL : (scalar *){fm};
#if TREE && _MPI
  restore_mpi (fp, slist);
#else // ! (TREE && _MPI)
#if !_MPI
  int rootlevel = 0;
#endif
#if TREE
  foreach_dimension()
    while ((1 << rootlevel) < header.n.x)
      rootlevel++;
  if (rootlevel > 0)
    init_grid (1 << rootlevel);
#endif // TREE
#if _MPI  
  foreach_cell() {
#else
  foreach_cell_restore (header.n, rootlevel) {
#endif
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    // skip subtree size
    fseek (fp, sizeof(double), SEEK_CUR);
    for (scalar s in slist) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
	fprintf (ferr, "restore(): error: expecting a scalar\n");
	exit (1);
      }
      if (s.i != INT_MAX)
	s[] = isfinite(val) ? val : nodata;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }
#if _GPU
  for (scalar s in slist)
    if (s.i != INT_MAX)
      s.gpu.stored = 1; // stored on CPU
#endif // _GPU
  for (scalar s in all)
    s.dirty = true;
#endif // ! (TREE && _MPI)
  
  scalar * other = NULL;
  for (scalar s in all)
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  free (other);
  
  free (slist);
  if (file)
    fclose (fp);

  // the events are advanced to catch up with the time  
  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);
  
  return true;
}

#endif // MULTIGRID

#if _GPU
# include "grid/gpu/output.h"
#endif
