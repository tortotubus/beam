#pragma once

#include <stdlib.h>
#include <string.h> /* memcpy */
#include <stddef.h> /* size_t */

/*!
 * \brief Immersed-boundary node state.
 *
 * Ownership:
 * - \c stencil.p is owned by the \c IBNode and must be freed with \ref
 * ibnode_free (or \ref ibnode_free_stencil).
 */
typedef struct {
  coord lagpos; /*!< Lagrangian nodal position */
  coord force;  /*!< Lagrange multiplier as a force */
  coord lagvel; /*!< Lagrangian nodal velocity */
  coord eulvel; /*!< Eulerian fluid velocity interpolated at the Lagrangian node */

  coord rhs; /*!< b_i = lagvel_i - (J u*)_i (constraint mismatch) */
  coord res; /*!< r_i (residual) */
  coord w;   /*!< w_i (search direction) */
  coord Ay;  /*!< y_i = (A w)_i */

  double weight; /*!< Quadrature weight */
  coord gravity; /*!< Gravity */

  int pid; /*!< MPI rank that owns this node */

  Cache stencil; /*!< Stencil cache (owns stencil.p) */
} IBNode;

/*!
 * \brief Compute the stencil storage length for the current Basilisk dimension.
 * \return Number of \c Index entries required for the stencil buffer.
 */
size_t ibnode_stencil_capacity_ (void) {
#if dimension == 1
  return (size_t) ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#elif dimension == 2
  return (size_t) sq ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#else
  return (size_t) cube ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#endif
}

/*!
 * \brief Initialize the \c stencil member of an \c IBNode.
 *
 * \param node Node to initialize.
 * \return 0 on success, -1 on allocation failure or invalid input.
 *
 * Allocates \c node->stencil.p and sets:
 * - \c stencil.n = 0
 * - \c stencil.nm = required capacity for the stencil buffer
 */
int ibnode_init_stencil (IBNode* node) {
  if (!node)
    return -1;

  node->stencil.n = 0;
  node->stencil.nm = ibnode_stencil_capacity_ ();
  node->stencil.p = NULL;

  if (node->stencil.nm == 0)
    return 0;

  node->stencil.p = (Index*) malloc (node->stencil.nm * sizeof (Index));
  if (!node->stencil.p) {
    node->stencil.n = 0;
    node->stencil.nm = 0;
    return -1;
  }
  return 0;
}

/*!
 * \brief Free the \c stencil member of an \c IBNode.
 *
 * \param node Node whose stencil should be freed.
 *
 * After this call:
 * - \c stencil.p == NULL
 * - \c stencil.n == 0
 * - \c stencil.nm == 0
 */
void ibnode_free_stencil (IBNode* node) {
  if (!node)
    return;
  free (node->stencil.p);
  node->stencil.p = NULL;
  node->stencil.n = 0;
  node->stencil.nm = 0;
}

/*!
 * \brief Initialize an \c IBNode.
 *
 * \param node Node to initialize.
 * \return 0 on success, -1 on allocation failure or invalid input.
 *
 * Sets all vector-like members to zero, sets default scalars, and initializes
 * the stencil buffer.
 */
int ibnode_init (IBNode* node) {
  if (!node)
    return -1;

  node->lagpos.x = 0.;
  node->force.x = 0.;
  node->lagvel.x = 0.;
  node->eulvel.x = 0.;
  node->rhs.x = 0.;
  node->res.x = 0.;
  node->w.x = 0.;
  node->Ay.x = 0.;
  node->gravity.x = 0.;

  node->lagpos.y = 0.;
  node->force.y = 0.;
  node->lagvel.y = 0.;
  node->eulvel.y = 0.;
  node->rhs.y = 0.;
  node->res.y = 0.;
  node->w.y = 0.;
  node->Ay.y = 0.;
  node->gravity.y = 0.;

  node->lagpos.z = 0.;
  node->force.z = 0.;
  node->lagvel.z = 0.;
  node->eulvel.z = 0.;
  node->rhs.z = 0.;
  node->res.z = 0.;
  node->w.z = 0.;
  node->Ay.z = 0.;
  node->gravity.z = 0.;

  node->weight = 1.;
  node->pid = -1;

  return ibnode_init_stencil (node);
}

/*!
 * \brief Free resources owned by an \c IBNode.
 *
 * \param node Node to free.
 *
 * Currently this frees only owned dynamic resources (the stencil buffer).
 */
void ibnode_free (IBNode* node) {
  if (!node)
    return;
  ibnode_free_stencil (node);
}

/*!
 * \brief Deep-copy an \c IBNode.
 *
 * \param dst Destination node (will be overwritten; any owned resources are
 * freed first).
 * \param src Source node.
 * \return 0 on success, -1 on allocation failure or invalid input.
 *
 * This performs a deep copy of \c stencil.p, so \c dst and \c src can be freed
 * independently.
 */
int ibnode_copy (IBNode* dst, const IBNode* src) {
  if (!dst || !src)
    return -1;

  /* Free any resources currently owned by dst */
  ibnode_free (dst);

  /* Copy trivially-copiable fields first (but DO NOT keep src->stencil.p) */
  *dst = *src;
  dst->stencil.p = NULL;

  if (src->stencil.nm > 0) {
    dst->stencil.p =
      (Index*) malloc ((size_t) src->stencil.nm * sizeof (Index));
    if (!dst->stencil.p) {
      dst->stencil.n = 0;
      dst->stencil.nm = 0;
      return -1;
    }
    memcpy (dst->stencil.p,
            src->stencil.p,
            (size_t) src->stencil.nm * sizeof (Index));
  }

  return 0;
}

/*!
 * \brief Move an \c IBNode (transfer ownership).
 *
 * \param dst Destination node (will be overwritten; any owned resources are
 * freed first).
 * \param src Source node (left in a valid, empty state).
 *
 * After the move, \c src does not own any stencil memory.
 */
void ibnode_move (IBNode* dst, IBNode* src) {
  if (!dst || !src)
    return;

  ibnode_free (dst);

  *dst = *src;

  src->stencil.p = NULL;
  src->stencil.n = 0;
  src->stencil.nm = 0;
}

/*!
 * \brief Swap two \c IBNode values.
 *
 * \param a First node.
 * \param b Second node.
 *
 * This is safe under the ownership model above (it swaps the owned pointer).
 */
void ibnode_swap (IBNode* a, IBNode* b) {
  if (!a || !b)
    return;
  IBNode tmp = *a;
  *a = *b;
  *b = tmp;
}
