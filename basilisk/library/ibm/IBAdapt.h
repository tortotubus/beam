#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBMacros.h"
#include "library/ibm/IBLocate.h"

struct Adapt2 {
  scalar* slist; // list of scalars
  double* max;   // tolerance for each scalar
  int* maxlevel; // maximum level of refinement for each scalar
  int minlevel;  // minimum level of refinement (default 1)
  scalar* list;  // list of fields to update (default all)
};

// scalar ib_noise_1[];
// scalar ib_noise_2[];

astats adapt_wavelet_ibm (scalar* slist,
                          double* max,
                          int maxlevel,
                          int minlevel = 1,
                          scalar* list = all);

astats adapt_wavelet2 (scalar* slist,
                       double* max,
                       int* maxlevel,
                       int minlevel = 1,
                       scalar* list = all);

scalar ib_noise_0[];

astats adapt_wavelet_ibm (scalar* slist,
                          double* max,
                          int maxlevel,
                          int minlevel = 1,
                          scalar* list = all) {
  int iblevel_0 = 0;
  // scalar ib_noise_1[]; int iblevel_1;
  // scalar ib_noise_2[]; int iblevel_2;

  foreach_ibnode_per_ibmesh () {
    if (mesh->depth > iblevel_0)
      iblevel_0 = mesh->depth;
  }

  scalar* slist_c = slist ? list_copy (slist) : NULL;
  slist_c = list_append (slist_c, ib_noise_0);
  int n = list_len (slist_c);

  double* max_c = (double*) malloc ((size_t) n * sizeof (double));
  int* maxlevel_c = (int*) malloc ((size_t) n * sizeof (int));
  assert (max_c && maxlevel_c);

  for (int i = 0; i < n - 1; i++) {
    max_c[i] = max[i];
    maxlevel_c[i] = maxlevel;
  }
  max_c[n - 1] = 1e-6;
  maxlevel_c[n - 1] = iblevel_0;

  foreach_cell () {
    ib_noise_0[] = 0.;
  }
  boundary ({ib_noise_0});

  foreach_ibnode_per_ibmesh () {
    int r = PESKIN_SUPPORT_RADIUS;
    double d = L0 / (1 << node->depth);
#if dimension == 1
    for (int i = -r; i <= r; i++) {
      coord e = {pos.x + i * d, pos.y, pos.z};
      coord_periodic_boundary (e);
      Point point = locate_nonlocal (e.x, e.y, e.z);
      int ig = 0, jg = 0, kg = 0;
      NOT_UNUSED (ig);
      NOT_UNUSED (jg);
      NOT_UNUSED (kg);
      if (point.level >= 0) {
        double rd = sq (x - c.x);
        ib_noise_0[] += exp (-rd / (2 * Delta));
      }
    }
#elif dimension == 2
    for (int i = -r; i <= r; i++) {
      for (int j = -r; j <= r; j++) {
        coord e = {pos.x + i * d, pos.y + j * d, pos.z};
        coord_periodic_boundary (e);
        Point point = locate_nonlocal (e.x, e.y, e.z);
        int ig = 0, jg = 0, kg = 0;
        NOT_UNUSED (ig);
        NOT_UNUSED (jg);
        NOT_UNUSED (kg);
        POINT_VARIABLES ();
        if (point.level >= 0) {
          double rd = sq (x - pos.x) + sq (y - pos.y);
          ib_noise_0[] += exp (-rd / (2 * Delta));
        }
      }
    }
#else
    for (int i = -r; i <= r; i++) {
      for (int j = -r; j <= r; j++) {
        for (int k = -r; k <= r; k++) {
          coord e = {pos.x + i * d, pos.y + j * d, pos.z + k * d};
          coord_periodic_boundary (e);
          Point point = locate_nonlocal (e.x, e.y, e.z);
          int ig = 0, jg = 0, kg = 0;
          NOT_UNUSED (ig);
          NOT_UNUSED (jg);
          NOT_UNUSED (kg);
          if (point.level >= 0) {
            double rd = sq (x - c.x) + sq (y - c.y) + sq (z - c.z);
            ib_noise_0[] += exp (-rd / (2 * Delta));
          }
        }
      }
    }
#endif
  }
  boundary ({ib_noise_0});

  astats st = adapt_wavelet2 (slist_c, max_c, maxlevel_c, minlevel, list);

  free (slist_c);
  free (max_c);
  free (maxlevel_c);

  ibmm.dirty = true;
  ibmeshmanager_update_pid ();

  return st;
}

trace astats adapt_wavelet2 (scalar* slist,
                             double* max,
                             int* maxlevel,
                             int minlevel = 1,
                             scalar* list = all) {
  scalar* ilist = list;

  if (is_constant (cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (slist);
  } else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
        list = list_add (list, s);
    }
    boundary (list);
    scalar* listr = list_concat ({cm}, slist);
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar* listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);

  // refinement
  if (minlevel < 1)
    minlevel = 1;

  int overall_maxlevel = minlevel;
  if (slist && maxlevel)
    for (int i = 0; i < list_len (slist); i++)
      if (maxlevel[i] > overall_maxlevel)
        overall_maxlevel = maxlevel[i];

  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell () {
    if (is_active (cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        }
        continue;
      } else { // !is_leaf (cell)
        if (cell.flags & refined) {
          // cell has already been refined, skip its children
          cell.flags &= ~too_coarse;
          continue;
        }
        // check whether the cell or any of its children is local
        bool local = is_local (cell);
        if (!local) {
          foreach_child () {
            if (is_local (cell)) {
              local = true;
              break;
            }
          }
        }
        if (local) {
          int i = 0;
          static const int just_fine = 1 << (user + 3);
          for (scalar s in slist) {
            double emax = max[i], sc[(1 << dimension) * s.block];
            int mlev = maxlevel[i++];
            double* b = sc;
            foreach_child () foreach_blockf (s)* b++ = s[];
            s.prolongation (point, s);
            b = sc;
            foreach_child () {
              foreach_blockf (s) {
                double e = fabs (*b - s[]);
                if (e > emax && level < mlev) {
                  cell.flags &= ~too_fine;
                  cell.flags |= too_coarse;
                } else if ((e <= emax / 1.5 || level > mlev) &&
                           !(cell.flags & (too_coarse | just_fine))) {
                  if (level >= minlevel)
                    cell.flags |= too_fine;
                } else if (!(cell.flags & too_coarse)) {
                  cell.flags &= ~too_fine;
                  cell.flags |= just_fine;
                }
                s[] = *b++;
              }
            }
          }
          foreach_child () {
            cell.flags &= ~just_fine;
            if (!is_leaf (cell)) {
              cell.flags &= ~too_coarse;
              if (level >= overall_maxlevel)
                cell.flags |= too_fine;
            } else if (!is_active (cell))
              cell.flags &= ~too_coarse;
          }
        }
      }
    } else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth (); l >= 0; l--) {
    foreach_cell () if (!is_boundary (cell)) {
      if (level == l) {
        if (!is_leaf (cell)) {
          if (cell.flags & refined)
            // cell was refined previously, unset the flag
            cell.flags &= ~(refined | too_fine);
          else if (cell.flags & too_fine) {
            if (is_local (cell) && coarsen_cell (point, listc))
              st.nc++;
            cell.flags &= ~too_fine; // do not coarsen parent
          }
        }
        if (cell.flags & too_fine)
          cell.flags &= ~too_fine;
        else if (level > 0 && (aparent (0).flags & too_fine))
          aparent (0).flags &= ~too_fine;
        continue;
      } else if (is_leaf (cell))
        continue;
    }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);
  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    free (list);

  return st;
}
