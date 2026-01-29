#pragma once

#include "library/ibm/kernels.h"
#include <stdlib.h>

/*!
 * @brief Initialize an IBNode struct
 */

typedef struct {
  coord lagpos; /*!< lagrangian nodal position */
  coord force;  /*!< lagrange multiplier as a force */
  coord lagvel; /*!< lagrangian nodal velocity */
  coord
    eulvel; /*!< eulerian fluid velocity interpolated at the lagrangian node */

  coord rhs; /*!< b_i = lagvel_i - (Ju*)_i (constraint mismatch) */
  coord res; /*!< r_i (residual) */
  coord w;   /*!< w_i (search direction) */
  coord Ay;  /*!< y_i = (A w)_i */

  double weight; /*!< quadtrature weight */
  coord gravity; /*!< gravity */

  int pid; /*!< MPI pid */

  Cache stencil; /*!< stencil cache */
} IBNode;

/*
 * Forward declaration
 */

void ibnode_init_stencil (IBNode* node);
void ibnode_free_stencil (IBNode* node);
void ibnode_init (IBNode* node);
void ibnode_free (IBNode* node);

/*
 * Actual definitions
 */

/*!
 * @brief Constructor for the stencil member
 *
 * @memberof IBNode
 */
void ibnode_init_stencil (IBNode* node) {
  node->stencil.n = 0;

#if dimension == 1
  size_t stencil_size = ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#elif dimension == 2
  size_t stencil_size = sq ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#else
  size_t stencil_size = cube ((2 * PESKIN_SUPPORT_RADIUS) + 1);
#endif

  node->stencil.nm = stencil_size;
  node->stencil.p = (Index*) malloc (stencil_size * sizeof (Index));
}

/*!
 * @brief Desctructor for the stencil member
 * @memberof IBNode
 */
void ibnode_free_stencil (IBNode* node) {
  free (node->stencil.p);
}

/*!
 * @brief Constructor for the IBNode struct
 * @memberof IBNode
 */
void ibnode_init (IBNode* node) {
  foreach_dimension () {
    node->lagpos.x = 0;
    node->force.x = 0;
    node->lagvel.x = 0;
    node->eulvel.x = 0;
    node->rhs.x = 0;
    node->res.x = 0;
    node->w.x = 0;
    node->Ay.x = 0;
    node->gravity.x = 0;
  }

  node->weight = 1.;
  node->pid = -1;

  ibnode_init_stencil (node);
}

/*!
 * @brief Destructor for the IBNode struct
 * @memberof IBNode
 */
void ibnode_free (IBNode* node) {
  ibnode_free_stencil (node);
}

void ibnode_copy (IBNode* to, IBNode* from) {
  *to = *from;
}

void ibnode_swap (IBNode* a, IBNode* b) {
  IBNode tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}
