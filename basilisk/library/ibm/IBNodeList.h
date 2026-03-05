#include "library/ibm/IBNode.h"

/* Type Definitions */

/**
 * @struct IBNodeList
 * @brief Non-owning dynamic array of pointers to @ref IBNode.
 *
 * @details The container stores borrowed pointers; it does not allocate or
 * destroy the referenced IBNode objects. It simply manages the pointer array.
 */
typedef struct {
  IBNode** ptrs;   /**< borrowed pointers to nodes */
  size_t size;     /**< active count */
  size_t capacity; /**< allocated length of ptrs[] */
} IBNodeList;

/* Function Declarations */

int ibnodelist_init (IBNodeList* list, size_t initial_capacity);
void ibnodelist_free (IBNodeList* list);
int ibnodelist_reserve (IBNodeList* list, size_t newcap);
int ibnodelist_push (IBNodeList* list, IBNode* node);
IBNode* ibnodelist_remove_swap (IBNodeList* list, size_t i);
static inline IBNode* ibnodelist_get (IBNodeList* list, size_t i);
void ibnodelist_clear (IBNodeList* list);
static size_t ibnodelist_next_cap_ (size_t cap);

/* Function Definitions */

/**
 * @brief
 */
static size_t ibnodelist_next_cap_ (size_t cap) {
  return cap ? cap * 2u : 8u;
}

/**
 * @brief Initialize an IBNodeList with an optional initial capacity.
 *
 * @param list Pointer to the IBNodeList to initialize. Must not be NULL.
 * @param initial_capacity Initial number of pointers to allocate. If 0, no
 * memory is allocated.
 *
 * @return 0 on success, -1 on failure (invalid input or memory allocation
 * failure).
 *
 * @details Sets the size to 0 and capacity to 0. If initial_capacity > 0,
 * allocates memory for that many IBNode pointers.
 *
 * @relates IBNodeList
 */
int ibnodelist_init (IBNodeList* list, size_t initial_capacity) {
  if (!list)
    return -1;
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;

  if (initial_capacity == 0)
    return 0;
  list->ptrs = (IBNode**) malloc (initial_capacity * sizeof (*list->ptrs));
  if (!list->ptrs)
    return -1;
  list->capacity = initial_capacity;
  return 0;
}

/**
 * @brief Free all allocated memory and reset the IBNodeList.
 *
 * @param list Pointer to the IBNodeList to free. If NULL, this function does
 * nothing.
 *
 * @details Deallocates the internal pointer array, resets size and capacity to
 * 0, and sets ptrs to NULL. Does not free the IBNode objects themselves.
 *
 * @relates IBNodeList
 */
void ibnodelist_free (IBNodeList* list) {
  if (!list)
    return;
  free (list->ptrs);
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;
}

/**
 * @brief Reserve capacity for at least newcap elements in the IBNodeList.
 *
 * @param list Pointer to the IBNodeList. Must not be NULL.
 * @param newcap New capacity to reserve.
 *
 * @return 0 on success or if newcap <= current capacity, -1 on failure (memory
 * allocation error).
 *
 * @details If newcap is less than or equal to the current capacity, no action
 * is taken. Otherwise, reallocates the internal pointer array to accommodate
 * newcap elements.
 *
 * @relates IBNodeList
 */
int ibnodelist_reserve (IBNodeList* list, size_t newcap) {
  if (!list)
    return -1;
  if (newcap <= list->capacity)
    return 0;
  IBNode** p = (IBNode**) realloc (list->ptrs, newcap * sizeof (*p));
  if (!p)
    return -1;
  list->ptrs = p;
  list->capacity = newcap;
  return 0;
}

/**
 * @brief Add an IBNode pointer to the end of the IBNodeList.
 *
 * @param list Pointer to the IBNodeList. Must not be NULL.
 * @param node Pointer to the IBNode to add. Must not be NULL.
 *
 * @return 0 on success, -1 on failure (invalid input or memory allocation
 * failure).
 *
 * @details If the list is at capacity, automatically reserves more capacity
 * (doubling the current capacity, or setting to 8 if empty) before adding the
 * node.
 *
 * @relates IBNodeList
 */
int ibnodelist_push (IBNodeList* list, IBNode* node) {
  if (!list || !node)
    return -1;
  if (list->size == list->capacity) {
    size_t newcap = ibnodelist_next_cap_ (list->capacity);
    if (ibnodelist_reserve (list, newcap) != 0)
      return -1;
  }
  list->ptrs[list->size++] = node;
  return 0;
}

/**
 * @brief Remove an element at index i using swap-with-last removal.
 *
 * @param list Pointer to the IBNodeList. Must not be NULL.
 * @param i Index of the element to remove. Must be < size.
 *
 * @return Pointer to the removed IBNode, or NULL if invalid input.
 *
 * @details Swaps the element at index i with the last element and decrements
 * size. This is efficient but does not preserve list order. The caller is
 * responsible for deciding how or if to destroy the returned IBNode object.
 *
 * @relates IBNodeList
 */
IBNode* ibnodelist_remove_swap (IBNodeList* list, size_t i) {
  if (!list || i >= list->size)
    return NULL;
  IBNode* dead = list->ptrs[i];
  list->ptrs[i] = list->ptrs[list->size - 1];
  list->size--;
  return dead; /* caller decides how/if to destroy it */
}

/**
 * @brief Retrieve the IBNode pointer at index i.
 *
 * @param list Pointer to the IBNodeList. If NULL, returns NULL.
 * @param i Index of the element to retrieve. Must be < size.
 *
 * @return Pointer to the IBNode at index i, or NULL if invalid input or out of
 * bounds.
 *
 * @relates IBNodeList
 */
IBNode* ibnodelist_get (IBNodeList* list, size_t i) {
  return (list && i < list->size) ? list->ptrs[i] : NULL;
}

/**
 * @brief Clear the list without freeing its storage.
 *
 * @param list Pointer to the IBNodeList. If NULL, this function does nothing.
 *
 * @details Sets size to 0 and leaves capacity/ptrs unchanged. Does not free
 * any IBNode objects (non-owning container).
 *
 * @relates IBNodeList
 */
void ibnodelist_clear (IBNodeList* list) {
  if (!list)
    return;
  list->size = 0;
}