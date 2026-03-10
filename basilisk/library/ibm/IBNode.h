#include "library/ibm/IBConfig.h"
#include "library/ibm/IBMacros.h"

/**
 * @struct IBNode
 *
 * @brief Immersed-boundary node state.
 */
typedef struct {
  coord pos;     /**< Lagrangian nodal position */
  int depth; /**< Desired depth of the cell for TREE */
  int pid; /**< MPI rank that owns this node: -1 if unknown */
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

  node->pos.x = 0.;
  node->pos.y = 0.;
  node->pos.z = 0.;
  node->depth = 0; 

#if _MPI
  node->pid = -1;
#else 
  node->pid = 0;
#endif

  return 0;
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

void ibnode_update_pid (IBNode* node) {

  coord centre_coord = node->pos;
  coord_periodic_boundary (centre_coord);
  Point centre_point = locate (centre_coord.x, centre_coord.y, centre_coord.z);

#if TREE
  if (centre_point.level > 0) {
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
    POINT_VARIABLES ();
    node->pid = cell.pid;
  }
#endif
}
