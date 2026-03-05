#include "library/ibm/IBConfig.h"
#include "library/ibm/IBMacros.h"

/**
 * @struct IBNode
 *
 * @brief Immersed-boundary node state.
 */
typedef struct {
  coord lagpos;  /**< Lagrangian nodal position */
  coord force;   /**< Lagrange multiplier as a force */
  coord lagvel;  /**< Lagrangian nodal velocity */
  coord eulvel;  /**< Eulerian fluid velocity interpolated at the Lagrangian node */
  double weight; /**< \f(\sum w^2\f)  where \f(w\f) is the dimensionless weight for a cell */
  int pid;       /**< MPI rank that owns this node: -1 if unknown */
#if TREE
  Cache stencil; /**< Stencil cache (owns stencil.p) */
#endif
} IBNode;

/**
 * @brief Compute the stencil storage length for the current Basilisk dimension.
 * @return Number of @c Index entries required for the stencil buffer.
 *
 * @relates IBNode
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

#if TREE
/**
 * @brief Initialize the stencil member of an @ref IBNode.
 *
 * @param node Node to initialize.
 * @return 0 on success, -1 on allocation failure or invalid input.
 *
 * Allocates @c node->stencil.p and sets:
 * - @c stencil.n = 0
 * - @c stencil.nm = required capacity for the stencil buffer
 *
 * @relates IBNode
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
#endif 

#if TREE
/**
 * @brief Free the \c stencil member of an \c IBNode.
 *
 * @param node Node whose stencil should be freed.
 *
 * @relates IBNode
 */
void ibnode_free_stencil (IBNode* node) {
  if (!node)
    return;
  free (node->stencil.p);
  node->stencil.p = NULL;
  node->stencil.n = 0;
  node->stencil.nm = 0;
}
#endif

/**
 * @brief Initialize an \c IBNode.
 *
 * @param node Node to initialize.
 * @return 0 on success, -1 on allocation failure or invalid input.
 *
 * @relates IBNode
 */
int ibnode_init (IBNode* node) {
  if (!node)
    return -1;

  node->lagpos.x = 0.;
  node->force.x = 0.;
  node->lagvel.x = 0.;
  node->eulvel.x = 0.;

  node->lagpos.y = 0.;
  node->force.y = 0.;
  node->lagvel.y = 0.;
  node->eulvel.y = 0.;

  node->lagpos.z = 0.;
  node->force.z = 0.;
  node->lagvel.z = 0.;
  node->eulvel.z = 0.;

  node->weight = 0.;
  node->pid = -1;

#if TREE
  return ibnode_init_stencil (node);
#else 
  return 0;
#endif
}

/**
 * @brief Free resources owned by an \c IBNode.
 *
 * @param node Node to free.
 *
 * @relates IBNode
 */
void ibnode_free (IBNode* node) {
  if (!node)
    return;
#if TREE
  ibnode_free_stencil (node);
#endif 
}

/**
 * @brief Deep-copy an @ref IBNode.
 *
 * @param dst Destination node (will be overwritten; any owned resources are
 * freed first).
 * @param src Source node.
 * @return 0 on success, -1 on allocation failure or invalid input.
 *
 * @relates IBNode
 */
int ibnode_copy (IBNode* dst, const IBNode* src) {
  if (!dst || !src)
    return -1;

  /* Free any resources currently owned by dst */
  ibnode_free (dst);

  /* Copy trivially-copiable fields first (but DO NOT keep src->stencil.p) */
  *dst = *src;

#if TREE
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
#endif 

  return 0;
}

/**
 * @brief Move an @ref IBNode (transfer ownership).
 *
 * @param dst Destination node (will be overwritten; any owned resources are
 * freed first).
 * @param src Source node (left in a valid, empty state).
 *
 * @relates IBNode
 */
void ibnode_move (IBNode* dst, IBNode* src) {
  if (!dst || !src)
    return;

  ibnode_free (dst);

  *dst = *src;

#if TREE
  src->stencil.p = NULL;
  src->stencil.n = 0;
  src->stencil.nm = 0;
#endif 
}

/**
 * @brief Swap two @ref IBNode values.
 *
 * @param a First node.
 * @param b Second node.
 *
 * @relates IBNode
 */
void ibnode_swap (IBNode* a, IBNode* b) {
  if (!a || !b)
    return;
  IBNode tmp = *a;
  *a = *b;
  *b = tmp;
}

/**
 * @brief Populates the @c stencil member of @ref IBNode with adjacent cells for
 * its stencil.
 */
// void ibnode_fill_stencil_cache (IBNode* node) {
//   foreach_neighborhood_coord_level(node->lagpos, PESKIN_SUPPORT_RADIUS, ) {
//     cache_append(&node->stencil, point, 0);
//   }
// }

void ibnode_update_pid(IBNode *node) {

  coord centre_coord = node->lagpos;

  Point centre_point = locate_nonlocal(centre_coord.x,centre_coord.y,centre_coord.z);
  
  if (centre_point.level > 0)
  {
    Point point = {0};
    point = centre_point;
#if dimension == 1
    int ig = point.i;
#elif dimension == 2
    int ig = point.i;
    int jg = point.j;
#else 
    int ig = point.i;
    int jg = point.j;
    int kg = point.j;
#endif
    POINT_VARIABLES();
    node->pid = cell.pid;
  }
}
