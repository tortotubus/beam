#include "library/ibm/IBNodeList.h"

/**
 * @struct IBExchangeList
 * @brief Non-owning list of nodes exchanged with one peer MPI rank.
 */
typedef struct {
  int pid;        /**< peer rank, -1 until first insertion */
  double * buffs;
  size_t nbuff;
  IBNodeList nodes;
} IBExchangeList;

/* Function declarations */
int ibexchangelist_init (IBExchangeList * list, size_t initial_capacity);
void ibexchangelist_free (IBExchangeList * list);
void ibexchangelist_clear (IBExchangeList * list);
int ibexchangelist_set_pid (IBExchangeList * list, int pid);
int ibexchangelist_push (IBExchangeList * list, IBNode * node);
int ibexchangelist_push_unique (IBExchangeList * list, IBNode * node);
void ibexchangelist_init_buffer(IBExchangeList * list, size_t nbuff);
void ibexchangelist_free_buffer(IBExchangeList * list);

/* Function definitions */

int
ibexchangelist_init (IBExchangeList * list, size_t initial_capacity) {
  if (!list)
    return -1;
  list->pid = -1;
  return ibnodelist_init (&list->nodes, initial_capacity);
}

void
ibexchangelist_free (IBExchangeList * list) {
  if (!list)
    return;
  ibnodelist_free (&list->nodes);
  list->pid = -1;
}

void
ibexchangelist_clear (IBExchangeList * list) {
  if (!list)
    return;
  ibnodelist_clear (&list->nodes);
  list->pid = -1;
}

int
ibexchangelist_set_pid (IBExchangeList * list, int pid) {
  if (!list)
    return -1;
  if (list->pid < 0) {
    list->pid = pid;
    return 0;
  }
  return (list->pid == pid) ? 0 : -1;
}

int
ibexchangelist_push (IBExchangeList * list, IBNode * node) {
  if (!list || !node)
    return -1;

  return ibnodelist_push (&list->nodes, node);
}

int
ibexchangelist_push_unique (IBExchangeList * list, IBNode * node) {
  if (!list || !node)
    return -1;

  /* Reverse scan keeps expected fast path when duplicates are recent. */
  for (size_t i = list->nodes.size; i > 0; i--) {
    if (list->nodes.ptrs[i - 1] == node)
      return 0;
  }

  return ibnodelist_push (&list->nodes, node);
}

void ibexchangelist_init_buffer(IBExchangeList * list, size_t nbuffs) {
  list->nbuff = nbuffs;
  int nn = list->nodes.size;
  list->buffs = (double*) calloc(nbuffs * nn, sizeof(double));
}

void ibexchangelist_free_buffer(IBExchangeList * list) {
  list->nbuff = 0;
  free(list->buffs);
}
