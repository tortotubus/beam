#pragma once

#include "library/ibm/IBNode.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

/*!
 * \brief Dynamic array of \c IBNode.
 *
 * Ownership / lifetime:
 * - The list owns the nodes stored in \c nodes[0..size-1].
 * - \ref ibnodelist_add_copy deep-copies a node into the list using \c ibnode_copy.
 * - \ref ibnodelist_free releases all nodes in the list using \c ibnode_free, then frees storage.
 *
 * Invariants:
 * - \c 0 <= size <= capacity
 * - If \c capacity == 0 then \c nodes may be NULL
 */
typedef struct {
  IBNode  *nodes;     /*!< Pointer to contiguous storage of length \c capacity. */
  size_t   size;      /*!< Number of valid/owned nodes in \c nodes[0..size-1]. */
  size_t   capacity;  /*!< Allocated length of \c nodes in elements. */
} IBNodeList;

/*!
 * \brief Compute the next capacity when growing a list.
 *
 * \param cap Current capacity.
 * \return Suggested next capacity (geometric growth).
 *
 * This uses a small non-zero base capacity to handle the \c cap==0 case.
 */
static size_t ibnodelist_next_cap_(size_t cap) {
  return (cap == 0) ? 8u : (cap * 2u);
}

/*!
 * \brief Initialize an \c IBNodeList.
 *
 * \param list Pointer to the list to initialize.
 * \param initial_capacity Initial storage capacity (may be 0).
 * \return 0 on success, -1 on allocation failure or invalid input.
 *
 * After successful initialization:
 * - \c list->size is 0
 * - \c list->capacity is \p initial_capacity
 * - \c list->nodes is allocated if \p initial_capacity > 0
 *
 * \note This function does not call \c ibnode_init on capacity slots; nodes are
 *       initialized on insertion.
 */
int ibnodelist_init(IBNodeList *list, size_t initial_capacity) {
  if (!list) return -1;

  list->nodes    = NULL;
  list->size     = 0;
  list->capacity = 0;

  if (initial_capacity == 0) return 0;

  list->nodes = (IBNode *)malloc(initial_capacity * sizeof(*list->nodes));
  if (!list->nodes) return -1;

  list->capacity = initial_capacity;
  return 0;
}

/*!
 * \brief Free all nodes and storage owned by a list.
 *
 * \param list Pointer to the list to free (may be NULL).
 *
 * For each element in \c nodes[0..size-1], calls \c ibnode_free, then frees the
 * underlying storage. Resets \c nodes to NULL and \c size/capacity to 0.
 */
void ibnodelist_free(IBNodeList *list) {
  if (!list) return;

  for (size_t i = 0; i < list->size; i++) {
    ibnode_free(&list->nodes[i]);
  }

  free(list->nodes);
  list->nodes    = NULL;
  list->size     = 0;
  list->capacity = 0;
}

/*!
 * \brief Ensure the list has at least the given capacity.
 *
 * \param list Pointer to the list.
 * \param new_capacity Minimum required capacity.
 * \return 0 on success, -1 on allocation failure or invalid input.
 *
 * If \p new_capacity is less than or equal to the current capacity, this is a
 * no-op. Otherwise the underlying storage is reallocated.
 *
 * \warning Reallocation may move the storage; any outstanding pointers into
 *          \c list->nodes become invalid after this call.
 */
int ibnodelist_reserve(IBNodeList *list, size_t new_capacity) {
  if (!list) return -1;
  if (new_capacity <= list->capacity) return 0;

  IBNode *p = (IBNode *)realloc(list->nodes, new_capacity * sizeof(*p));
  if (!p) return -1;

  list->nodes    = p;
  list->capacity = new_capacity;
  return 0;
}

/*!
 * \brief Append a node to the list by deep-copying it.
 *
 * \param list Pointer to the list.
 * \param node Pointer to the node to copy from.
 * \return The inserted index on success, -1 on failure.
 *
 * This function grows the list if necessary, initializes the destination slot
 * via \c ibnode_init, then deep-copies the contents using \c ibnode_copy.
 *
 * Ownership:
 * - The list owns the copied node.
 * - The caller retains ownership of \p node and may free/modify it independently.
 *
 * \note This assumes \c ibnode_copy performs a deep copy compatible with
 *       \c ibnode_free.
 */
int ibnodelist_add_copy(IBNodeList *list, const IBNode *node) {
  if (!list || !node) return -1;

  if (list->size == list->capacity) {
    size_t new_cap = ibnodelist_next_cap_(list->capacity);
    if (new_cap < list->size + 1) new_cap = list->size + 1;
    if (ibnodelist_reserve(list, new_cap) != 0) return -1;
  }

  size_t idx = list->size++;
  ibnode_init(&list->nodes[idx]);        /* destination in known state */
  ibnode_copy(&list->nodes[idx], node);  /* deep copy into owned storage */

  return (idx <= (size_t)INT32_MAX) ? (int)idx : -1;
}

/*!
 * \brief Get a mutable pointer to an element by index.
 *
 * \param list Pointer to the list.
 * \param index Index in \c [0, size).
 * \return Pointer to the element, or NULL if out of range or \p list is NULL.
 */
IBNode *ibnodelist_get(IBNodeList *list, size_t index) {
  if (!list) return NULL;
  if (index >= list->size) return NULL;
  return &list->nodes[index];
}

/*!
 * \brief Get a const pointer to an element by index.
 *
 * \param list Pointer to the list.
 * \param index Index in \c [0, size).
 * \return Const pointer to the element, or NULL if out of range or \p list is NULL.
 */
const IBNode *ibnodelist_get_const(const IBNodeList *list, size_t index) {
  if (!list) return NULL;
  if (index >= list->size) return NULL;
  return &list->nodes[index];
}
