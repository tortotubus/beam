
#include "grid/events.h"

typedef struct {
  int iter, inext;
  double t, tnext;
  int n_ev;
  Event* events;
} EventState;

EventState get_event_state (void) {
  EventState es = {.iter = iter,
                   .inext = inext,
                   .t = t,
                   .tnext = tnext,
                   .n_ev = 0,
                   .events = NULL};

  // count events
  int n_ev = 0;
  for (Event* ev = Events; !ev->last; ev++)
    n_ev++;
  es.n_ev = n_ev;

  // allocate and copy stateful fields
  es.events = (Event*) calloc (n_ev, sizeof (Event));
  assert (es.events);

  int i_ev = 0;
  for (Event* ev = Events; !ev->last; ev++, i_ev++) {
    // Copy only the minimal stateful parts
    es.events[i_ev].i = ev->i;
    es.events[i_ev].t = ev->t;
    es.events[i_ev].a = ev->a;

    // Copy name safely (shallow is fine since restored by order)
    if (ev->name)
      es.events[i_ev].name = strdup (ev->name);
    else
      es.events[i_ev].name = NULL;
  }

  return es;
}

/**
 * Restore the global scheduler and per-event state from a saved snapshot.
 */
void set_event_state (const EventState* es) {
  assert (es);

  while (iter < es->iter && events (false))
    iter = inext;
  events (false);
  while (t < es->t && events (false))
    t = tnext;
  t = es->t;
  events (false);
  // t = es->t;
  // tnext = es->tnext;

  // // restore per-event state in order
  // int i_ev = 0;
  // for (Event *ev = Events; !ev->last && i_ev < es->n_ev; ev++, i_ev++) {
  //   ev->i = es->events[i_ev].i;
  //   ev->t = es->events[i_ev].t;
  //   ev->a = es->events[i_ev].a;
  // }

  // // recompute next event triggers
  // inext = END_EVENT;
  // tnext = HUGE;
  // for (Event *ev = Events; !ev->last; ev++) {
  //   if (ev->t > t && ev->t < tnext)
  //     tnext = ev->t;
  //   if (ev->i > iter && ev->i < inext)
  //     inext = ev->i;
  // }
}

#include "library/grid/morton.h"
#include "library/grid/foreach_cell_bfs.h"
#include "library/io/hdf/hdf.h"
#include "library/io/output-common.h"

typedef struct {
  double l0;
  double x0;
  double y0;
  double z0;
  int dim;
  int num_of_cells;
  int depth;
  double time;
  int iter;
  uint64_t* descriptors;
  unsigned short* flags;
} TreeTopology;

/**
 *
 *
 *
 *
 */
TreeTopology get_tree_topology () {
  TreeTopology tp;
  tp.l0 = L0;
  tp.x0 = X0;
  tp.y0 = Y0;
  tp.z0 = Z0;
  tp.dim = dimension;
  tp.num_of_cells = 0;
  tp.depth = grid->depth;
  tp.time = t;
  tp.iter = iter;

  foreach_cell_BFS () {
    tp.num_of_cells++;
  }

  tp.descriptors = qcalloc (tp.num_of_cells, uint64_t);
  tp.flags = qcalloc (tp.num_of_cells, unsigned short);

  uint64_t cell_iter = 0;
  foreach_cell_BFS () {
    tp.descriptors[cell_iter] = point_to_morton (point);
    tp.flags[cell_iter] = cell.flags;
    cell_iter++;
  }

  return tp;
}

/**
 *
 *
 *
 *
 */
void set_tree_topology (TreeTopology* tp) {
  init_grid (1 << 0);
  origin (tp->x0, tp->y0, tp->z0);
  size (tp->l0);
  scalar* slist = all;

  for (uint64_t di = 0; di < tp->num_of_cells; di++) {
    Point point = morton_to_point (tp->descriptors[di]);
    if (!((tp->flags[di]) & leaf))
      refine_cell (point, slist, 0, NULL);
  }

  update_cache ();
}

/**
 *
 *
 *
 *
 */
void dump_hdf (const char* basename, scalar* slist = all) {
  char fname[128];
  sprintf (fname, "%s", basename);

  // printf ("Writing checkpoint to %s\n", fname);
#if _MPI
  HDF hdf = HDF_init_MPIIO (fname, true);
#else
  HDF hdf = HDF_init (fname, true);
#endif

  /*
   * Create groups
   */

  hid_t grp_basilisk_id = H5Gcreate2 (
    hdf.file_id, "/basilisk", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t grp_grid_id = H5Gcreate2 (
    hdf.file_id, "/basilisk/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t grp_scalars_id = H5Gcreate2 (grp_grid_id,
                                     "/basilisk/grid/scalars",
                                     H5P_DEFAULT,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT);

  hid_t grp_events_id = H5Gcreate2 (
    grp_grid_id, "/basilisk/events", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /*
   * Get data objects
   */

  TreeTopology tp = get_tree_topology ();
  EventState es = get_event_state ();

  /*
   * Write /basilisk/grid attributes
   */
  HDF_write_scalar_attribute ("X0", &tp.x0, H5T_IEEE_F64LE, grp_grid_id, &hdf);
  HDF_write_scalar_attribute ("Y0", &tp.y0, H5T_IEEE_F64LE, grp_grid_id, &hdf);
  HDF_write_scalar_attribute ("Z0", &tp.z0, H5T_IEEE_F64LE, grp_grid_id, &hdf);
  HDF_write_scalar_attribute ("L0", &tp.l0, H5T_IEEE_F64LE, grp_grid_id, &hdf);
  HDF_write_scalar_attribute (
    "Dimensions", &tp.dim, H5T_NATIVE_INT, grp_grid_id, &hdf);
  HDF_write_scalar_attribute (
    "CellCount", &tp.num_of_cells, H5T_NATIVE_INT, grp_grid_id, &hdf);
  HDF_write_scalar_attribute (
    "Depth", &tp.depth, H5T_NATIVE_INT, grp_grid_id, &hdf);

  /*
   * Write /basilisk/events attributes
   */
  HDF_write_scalar_attribute (
    "iter", &es.iter, H5T_NATIVE_INT, grp_events_id, &hdf);
  HDF_write_scalar_attribute (
    "inext", &es.inext, H5T_NATIVE_INT, grp_events_id, &hdf);
  HDF_write_scalar_attribute ("t", &es.t, H5T_IEEE_F64LE, grp_events_id, &hdf);
  HDF_write_scalar_attribute (
    "tnext", &es.tnext, H5T_IEEE_F64LE, grp_events_id, &hdf);

  /*
   * Write /basilisk/events per-event arrays
   */
  {
    // pack arrays
    int* events_i = malloc (es.n_ev * sizeof (int));
    double* events_t = malloc (es.n_ev * sizeof (double));
    int* events_a = malloc (es.n_ev * sizeof (int));

    for (int k = 0; k < es.n_ev; k++) {
      events_i[k] = es.events[k].i;
      events_t[k] = es.events[k].t;
      events_a[k] = es.events[k].a;
    }

    // dataset dimensions
    hsize_t dims[1] = {es.n_ev};
    hsize_t chunk_dims[1] = {es.n_ev};

    // write compressed datasets
    HDF_write_compressed_dataset ("events_i",
                                  events_i,
                                  H5T_NATIVE_INT,
                                  grp_events_id,
                                  1,
                                  dims,
                                  dims,
                                  chunk_dims,
                                  0,
                                  &hdf);

    HDF_write_compressed_dataset ("events_t",
                                  events_t,
                                  H5T_IEEE_F64LE,
                                  grp_events_id,
                                  1,
                                  dims,
                                  dims,
                                  chunk_dims,
                                  0,
                                  &hdf);

    HDF_write_compressed_dataset ("events_a",
                                  events_a,
                                  H5T_NATIVE_INT,
                                  grp_events_id,
                                  1,
                                  dims,
                                  dims,
                                  chunk_dims,
                                  0,
                                  &hdf);

    // cleanup
    free (events_i);
    free (events_t);
    free (events_a);
  }

#if !_MPI

  // Descriptors
  {
    const char* descriptors_name = "Descriptors";
    uint64_t* descriptors_data = tp.descriptors;
    hid_t descriptors_dtype = H5T_STD_U64LE;
    hid_t descriptors_group = grp_grid_id;
    int descriptors_rank = 1;
    hsize_t descriptors_dims[] = {tp.num_of_cells};
    hsize_t descriptors_chunk_dims[] = {tp.num_of_cells};
    unsigned int descriptors_compression_level = 5;

    HDF_write_compressed_dataset (descriptors_name,
                                  descriptors_data,
                                  descriptors_dtype,
                                  descriptors_group,
                                  descriptors_rank,
                                  descriptors_dims,
                                  descriptors_dims,
                                  descriptors_chunk_dims,
                                  descriptors_compression_level,
                                  &hdf);
  }

  // Flags
  {
    const char* flags_name = "Flags";
    unsigned short* flags_data = tp.flags;
    hid_t flags_dtype = H5T_NATIVE_USHORT;
    hid_t flags_group = grp_grid_id;
    int flags_rank = 1;
    hsize_t flags_dims[] = {tp.num_of_cells};
    hsize_t flags_chunk_dims[] = {tp.num_of_cells};
    unsigned int flags_compression_level = 5;

    HDF_write_compressed_dataset (flags_name,
                                  flags_data,
                                  flags_dtype,
                                  flags_group,
                                  flags_rank,
                                  flags_dims,
                                  flags_dims,
                                  flags_chunk_dims,
                                  flags_compression_level,
                                  &hdf);
  }

  // Scalars
  for (scalar s in slist) {

    hsize_t scalar_chunk_dims[] = {tp.num_of_cells};
    hsize_t scalar_dims[] = {tp.num_of_cells};
    int scalar_rank = 1;

    // Value for the Scalar dataset
    double* s_data = malloc (tp.num_of_cells * sizeof (double));

    // Malloc error handling
    if (!s_data) {
      perror ("malloc(s_data)");
      exit (1);
    }

    // Copy the data from the tree into s_data
    for (int di = 0; di < tp.num_of_cells; di++) {
      Point point = morton_to_point (tp.descriptors[di]);
      s_data[di] = val (s);
    }

    HDF_write_compressed_dataset (s.name,            /* dataset_name */
                                  s_data,            /* data */
                                  H5T_IEEE_F64LE,    /* dtype_id */
                                  grp_scalars_id,    /* group_id */
                                  scalar_rank,       /* rank */
                                  scalar_dims,       /* dims[] */
                                  scalar_dims,       /* max_dims[] */
                                  scalar_chunk_dims, /* chunk_dims[] */
                                  5,                 /* compression_level */
                                  &hdf);
    free (s_data);

  } /* end of “for (scalar s in scalar_list)” */

  //

#endif

  H5Gclose (grp_grid_id);
  H5Gclose (grp_scalars_id);
  H5Gclose (grp_events_id);
  H5Gclose (grp_basilisk_id);

  HDF_close (&hdf);
}

/**
 *
 *
 *
 *
 */
void restore_hdf (const char* basename, scalar* slist = all) {
  char fname[128];
  sprintf (fname, "%s", basename);
#if !_MPI
  HDF hdf = HDF_open (fname);

  hid_t grp_basilisk_id = H5Gopen (hdf.file_id, "/basilisk", H5P_DEFAULT);
  hid_t grp_events_id =
    H5Gopen (grp_basilisk_id, "/basilisk/events", H5P_DEFAULT);
  hid_t grp_grid_id = H5Gopen (grp_basilisk_id, "/basilisk/grid", H5P_DEFAULT);
  hid_t grp_scalars_id =
    H5Gopen (grp_grid_id, "/basilisk/grid/scalars", H5P_DEFAULT);

  TreeTopology tp = {0};
  EventState es = {0};

  HDF_read_scalar_attribute ("X0", &tp.x0, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("Y0", &tp.y0, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("Z0", &tp.z0, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("L0", &tp.l0, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("Dimensions", &tp.dim, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("CellCount", &tp.num_of_cells, grp_grid_id, &hdf);
  HDF_read_scalar_attribute ("Depth", &tp.depth, grp_grid_id, &hdf);

  HDF_read_scalar_attribute ("iter", &es.iter, grp_events_id, &hdf);
  HDF_read_scalar_attribute ("inext", &es.inext, grp_events_id, &hdf);
  HDF_read_scalar_attribute ("t", &es.t, grp_events_id, &hdf);
  HDF_read_scalar_attribute ("tnext", &es.tnext, grp_events_id, &hdf);
  /*
   * Read /basilisk/events per-event arrays
   */
  {
    int* events_i = NULL;
    double* events_t = NULL;
    int* events_a = NULL;

    hsize_t dims[1] = {0};
    int rank = 1;

    // read datasets
    HDF_read_dataset ("events_i",
                      (void**) &events_i,
                      H5T_NATIVE_INT,
                      grp_events_id,
                      rank,
                      dims,
                      dims,
                      &hdf);

    HDF_read_dataset ("events_t",
                      (void**) &events_t,
                      H5T_IEEE_F64LE,
                      grp_events_id,
                      rank,
                      dims,
                      dims,
                      &hdf);

    HDF_read_dataset ("events_a",
                      (void**) &events_a,
                      H5T_NATIVE_INT,
                      grp_events_id,
                      rank,
                      dims,
                      dims,
                      &hdf);

    // store in EventState
    es.n_ev = (int) dims[0];
    es.events = calloc (es.n_ev, sizeof (Event));
    assert (es.events);

    for (int k = 0; k < es.n_ev; k++) {
      es.events[k].i = events_i[k];
      es.events[k].t = events_t[k];
      es.events[k].a = events_a[k];
    }

    free (events_i);
    free (events_t);
    free (events_a);
  }

  {
    const char* descriptors_name = "Descriptors";
    uint64_t* descriptors_data = NULL;
    hid_t descriptors_dtype = H5T_STD_U64LE;
    hid_t descriptors_group = grp_grid_id;
    int descriptors_rank = 1;
    hsize_t descriptors_dims[] = {0};

    HDF_read_dataset (descriptors_name,
                      (void**) &descriptors_data,
                      descriptors_dtype,
                      descriptors_group,
                      descriptors_rank,
                      descriptors_dims,
                      descriptors_dims,
                      &hdf);

    tp.descriptors = descriptors_data;
  }

  {
    const char* flags_name = "Flags";
    unsigned short* flags_data = NULL;
    hid_t flags_dtype = H5T_NATIVE_USHORT;
    hid_t flags_group = grp_grid_id;
    int flags_rank = 1;
    hsize_t flags_dims[] = {0};

    HDF_read_dataset (flags_name,
                      (void**) &flags_data,
                      flags_dtype,
                      flags_group,
                      flags_rank,
                      flags_dims,
                      flags_dims,
                      &hdf);

    tp.flags = flags_data;
  }

  set_tree_topology (&tp);

  for (scalar s in slist) {

    hsize_t scalar_dims[] = {tp.num_of_cells};
    int scalar_rank = 1;

    // Value for the Scalar dataset
    double* s_data = NULL;

    HDF_read_dataset (s.name,           /* dataset_name */
                      (void**) &s_data, /* data */
                      H5T_IEEE_F64LE,   /* dtype_id */
                      grp_scalars_id,   /* group_id */
                      scalar_rank,      /* rank */
                      scalar_dims,      /* dims[] */
                      scalar_dims,      /* max_dims[] */
                      &hdf);

    // Copy the data from the s_data into the tree
    for (int mi = 0; mi < tp.num_of_cells; mi++) {
      Point point = morton_to_point (tp.descriptors[mi]);
      val (s) = s_data[mi];
    }

    free (s_data);
  } /* end of “for (scalar s in scalar_list)” */

#endif

  set_event_state (&es);

  H5Gclose (grp_scalars_id);
  H5Gclose (grp_grid_id);
  H5Gclose (grp_events_id);
  H5Gclose (grp_basilisk_id);
  HDF_close (&hdf);
}
