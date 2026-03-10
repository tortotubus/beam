#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBMacros.h"
#include "library/ibm/IBLocate.h"
 

astats adapt_wavelet_ibm (scalar* slist,      // list of scalars
                          double* max,        // tolerance for each scalar
                          int maxlevel,       // maximum level of refinement
                          int minlevel = 1,   // minimum level of refinement
                          scalar* list = all) // list of fields to update
{
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
    scalar* listr = list_concat (slist, {cm});
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
  tree->refined.n = 0;

  static const int refined = 1 << user;
  static const int too_fine = 1 << (user + 1);
  static const int too_coarse = 1 << (user + 2);
  static const int just_fine = 1 << (user + 3);
  static const int just_fine_ibm = 1 << (user + 4);

  // ibm refinement
  for (int it = 0; it < 20; it++) {
    int nrf = 0, nrc = 0;
    foreach_ibnode() {
      coord c = node->pos;
      coord_periodic_boundary (c);
      int L = node->depth;
      int l = locate_nonlocal (c.x, c.y, c.z).level;
      if (l >= 0) {
        int d = L - l;
        // int r = (PESKIN_SUPPORT_RADIUS + (1 << d) - 1) >> d;
        int r;
        if (d > 0) {                                       // Too coarse
          r = (PESKIN_SUPPORT_RADIUS + (1 << d) - 1) >> d; // ceil((R+1)/2^d)
          // foreach_neighbor_coord_nonlocal (r, c) {
          foreach_neighborhood_coord_nonlocal(c,r) {
            if (is_local (cell)) {
              refine_cell (point, listc, refined, &tree->refined);
              st.nf++;
              nrf++;
            }
          }

        } else if (d < 0) { // Too fine
          r = (PESKIN_SUPPORT_RADIUS) << (-d);
          // foreach_neighbor_coord_nonlocal (r, c) {
          foreach_neighborhood_coord_nonlocal(c,r) {
            if (is_local (cell)) {
              coarsen_cell (point, listc);
              st.nc++;
              nrc++;
            }
          }

        } else { // Juuust right
          r = PESKIN_SUPPORT_RADIUS;
          // foreach_neighbor_coord_nonlocal (r, c) {
          foreach_neighborhood_coord_nonlocal(c,r) {
            if (is_local (cell)) {
              if ((aparent (0).flags & leaf)) {
                point.level--;
#if dimension == 1
                point.i = ((point.i - GHOSTS) >> 1) + GHOSTS;
#elif dimension == 2
                point.i = ((point.i - GHOSTS) >> 1) + GHOSTS;
                point.j = ((point.j - GHOSTS) >> 1) + GHOSTS;
#else
                point.i = ((point.i - GHOSTS) >> 1) + GHOSTS;
                point.j = ((point.j - GHOSTS) >> 1) + GHOSTS;
                point.k = ((point.k - GHOSTS) >> 1) + GHOSTS;
#endif
                refine_cell (point, listc, refined, &tree->refined);
                st.nf++;
                nrf++;
              } else {
                cell.flags |= just_fine_ibm;
              }

            }
          }
        }
      }
      // mpi_boundary_refine (listc);
    }
    mpi_boundary_refine (listc);
    mpi_boundary_update (listc);
    mpi_all_reduce (nrf, MPI_INT, MPI_SUM);
    mpi_all_reduce (nrc, MPI_INT, MPI_SUM);
    if (!nrf && !nrc)
      break; // reached target everywhere
  }

  // recurse
  foreach_cell () {
    if (is_active (cell)) {
      // Actual refinement
      if (is_leaf (cell)) {
        if ((cell.flags & too_coarse) && !(cell.flags & just_fine_ibm)) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        } else if ((cell.flags & too_coarse) && (cell.flags & just_fine_ibm)) {
          cell.flags &= ~too_coarse;
        }
        continue;

      } else { // !is_leaf (cell)

        // Tagging for refinement

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

          for (scalar s in slist) {
            double emax = max[i++], sc[1 << dimension];
            int c = 0;
            foreach_child () {
              sc[c++] = s[];
            }
            s.prolongation (point, s);
            c = 0;

            foreach_child () {
              double e = fabs (sc[c] - s[]);
              if (e > emax &&
                  (level < maxlevel || cell.flags & just_fine_ibm)) {
                cell.flags &= ~too_fine;
                cell.flags |= too_coarse;
              } else if ((e <= emax / 1.5 || level > maxlevel) &&
                         !(cell.flags & (too_coarse | just_fine))) {
                if (level >= minlevel)
                  cell.flags |= too_fine;
              } else if (!(cell.flags & too_coarse)) {
                cell.flags &= ~too_fine;
                cell.flags |= just_fine;
              }
              s[] = sc[c++];
            }
          }

          foreach_child () {
            cell.flags &= ~just_fine;

            if (!is_leaf (cell)) {
              cell.flags &= ~too_coarse;
              if (level >= maxlevel && !(cell.flags & just_fine_ibm)) {
                cell.flags |= too_fine;
              }
            } else if (!is_active (cell)) {
              cell.flags &= ~too_coarse;
            }
          }
        }
      }
    } else { // inactive cell
      continue;
    }
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
          else if ((cell.flags & too_fine)) {
            bool just_fine_ibm_children = false;
            foreach_child () {
              if (cell.flags & just_fine_ibm) {
                just_fine_ibm_children = true;
              }
            }
            if (is_local (cell) && !just_fine_ibm_children &&
                coarsen_cell (point, listc))
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

  ibmm.dirty = true;

  return st;
}