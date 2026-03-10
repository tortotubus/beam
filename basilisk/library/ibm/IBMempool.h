#include "grid/mempool.h"
#include "library/ibm/IBFields.h"
#include "library/ibm/IBNodeList.h"

/**
 * @struct IBMempool
 * 
 * @brief 
 */
typedef struct {
  Mempool* pool;   /**< Memory pool pointer */
  IBNodeList active; /**< list of pointers to active IBNodes */
  int len;           /**< The (1D) size of the array */
  size_t datasize;   /**< Extra bytes stored after each IBNode */
} IBMempool;


/* Function declarations */

IBMempool ibmempool_init (size_t pool_bytes, size_t datasize);
void ibmempool_free (IBMempool* ibmp);
IBNode* ibmempool_alloc_node (IBMempool* ibmp);
long ibmempool_index_of (IBMempool* ibmp, IBNode* node);
int ibmempool_free_node_ptr (IBMempool* ibmp, IBNode* node);
void ibmempool_free_node (IBMempool* ibmp, long i);
inline size_t round_up_multiple(size_t n, size_t a);
static inline size_t ibmempool_stride (const IBMempool* ibmp);


/* Macro definitions */

#define alignof(T) offsetof(struct { char c; T x; }, x)
/* Function definitions */

/**
 * @brief 
 * 
 * @relates IBMempool
 */
inline size_t round_up_multiple(size_t n, size_t a)
{
  return a ? ((n + a - 1) / a) * a : n;
}

/**
 * @brief 
 * 
 * @relates IBMempool
 */
static inline size_t ibmempool_stride (const IBMempool* ibmp)
{
  assert (ibmp);
  return round_up_multiple (sizeof (IBNode) + ibmp->datasize, 8);
}

/**
 * @brief Initialize an IBMempool with a specified pool size and initial capacity.
 *
 * @param pool_bytes Total size in bytes for the memory pool.
 *
 * @return An initialized IBMempool structure.
 *
 * @details Creates a new memory pool for IBNode objects and allocates an array to
 * track active nodes. The pool block size is rounded up to the system's max alignment.
 * The active list starts empty and grows dynamically as nodes are added.
 * 
 * @relates IBMempool
 */
IBMempool ibmempool_init (size_t pool_bytes, size_t datasize) {
  IBMempool ibmp = {0};

  ibmp.datasize = datasize;
  size_t block = ibmempool_stride (&ibmp);
  ibmp.pool = mempool_new (pool_bytes, block);

  ibnodelist_init (&ibmp.active, 0);

  return ibmp;
}

/**
 * @brief Destroy an IBMempool and free all associated resources.
 *
 * @param mp Pointer to the IBMempool to destroy. Must not be NULL.
 *
 * @details Deallocates the active nodes list and destroys the underlying memory pool.
 * Resets all pointers to NULL and counts to 0. Does not free individual IBNode objects
 * (they are managed by the memory pool).
 * 
 * @relates IBMempool
 */
void ibmempool_free (IBMempool* ibmp) {
  assert(ibmp);
  ibnodelist_free (&ibmp->active);
  mempool_destroy (ibmp->pool);
  ibmp->pool = NULL;
}

/**
 * @brief Create a new IBNode from the memory pool and add it to the active list.
 *
 * @param ibmp Pointer to the IBMempool. Must not be NULL.
 *
 * @return Pointer to the newly created IBNode.
 *
 * @details Allocates a new IBNode from the pool. The new node is added to the active
 * list. Program aborts if allocation fails.
 * 
 * @relates IBMempool
 */
IBNode* ibmempool_alloc_node (IBMempool* ibmp) {
  IBNode* node = (IBNode*) mempool_alloc0 (ibmp->pool);
  assert (node);

  ibnodelist_push (&ibmp->active, node);

  return node;
}

/**
 * @brief Find the index of an active IBNode pointer in the mempool.
 *
 * @param ibmp Pointer to the IBMempool. Must not be NULL.
 * @param node Pointer to the IBNode to locate. Must not be NULL.
 *
 * @return Index of node in active array, or -1 if not found.
 * 
 * @relates IBMempool
 */
long ibmempool_index_of (IBMempool* ibmp, IBNode* node) {
  if (!ibmp || !node)
    return -1;
  for (size_t i = 0; i < ibmp->active.size; i++) {
    if (ibmp->active.ptrs[i] == node)
      return i;
  }
  return -1;
}

/**
 * @brief Release and destroy an IBNode by pointer.
 *
 * @param ibmp Pointer to the IBMempool. Must not be NULL.
 * @param node Pointer to the IBNode to release. Must not be NULL.
 *
 * @return 0 on success, -1 if the node is not found.
 * 
 * @relates IBMempool
 */
int ibmempool_free_node_ptr (IBMempool* ibmp, IBNode* node) {
  long idx = ibmempool_index_of (ibmp, node);
  if (idx < 0)
    return -1;
  ibmempool_free_node (ibmp, idx);
  return 0;
}

/**
 * @brief Release and destroy an IBNode at the given index.
 *
 * @param ibmp Pointer to the IBMempool. Must not be NULL.
 * @param i Index of the node to release in the active array. Must be >= 0 and < size.
 *
 * @details Removes the node at index i from the active list using swap-with-last.
 * Calls ibnode_free() to clean up the node, then returns the memory to the pool
 * via mempool_free(). Index out of bounds causes program abort.
 * 
 * @relates IBMempool
 */
void ibmempool_free_node (IBMempool* ibmp, long i) {
  assert (i >= 0 && (size_t) i < ibmp->active.size);

  IBNode* dead = ibnodelist_remove_swap (&ibmp->active, (size_t) i);

  ibnode_free (dead);

  mempool_free (ibmp->pool, dead);
}
