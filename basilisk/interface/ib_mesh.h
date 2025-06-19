#include <stdlib.h>

/**
 * Distance functions
 */

#define ACROSS_PERIODIC(a,b) (fabs(a - b) > L0/2.)
#define PERIODIC_1DIST(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DIST(a,b) : a - b)
#define PERIODIC_1DAVG(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DAVG(a,b) : a + b)
#define GENERAL_SQNORM(a,b) (sq(GENERAL_1DIST(a.x, b.x)) + \
  sq(GENERAL_1DIST(a.y, b.y)) + sq(GENERAL_1DIST(a.z, b.z)))

/**
 * IBNode
 */
#if dimension == 1
  #define STENCIL_SIZE 5
#elif dimension == 2
  #define STENCIL_SIZE 25
#else // dimension == 3
  #define STENCIL_SIZE 125
#endif

typedef struct
{
  coord pos;
  coord velocity;
  coord force;
  Cache stencil;
} IBNode;

void init_ibnode_stencil(IBNode *n) {
  n->stencil.n = STENCIL_SIZE;
  n->stencil.nm = STENCIL_SIZE;
  n->stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));  
}
void free_ibnode_stencil(IBNode *n) {
  free(n->stencil.p);
}

void init_ibnode(IBNode *n) {
  init_ibnode_stencil(n);
}

void free_ibnode(IBNode *n) {
  free_ibnode_stencil(n);
}

/**
 * IBMesh
 */

typedef struct
{
  IBNode* nodes;
  size_t n_nodes;
} IBMesh;

void init_ibmesh(IBMesh* m, size_t n_nodes) {
  m->nodes = (IBNode*) calloc(n_nodes, sizeof(IBNode));
  m->n_nodes = n_nodes;
  for (size_t i = 0; i < n_nodes; i++) {
    init_ibnode(&m->nodes[i]);
  }
}

void free_ibmesh(IBMesh* m)
{
  if (m->nodes) {
    for (size_t i = 0; i < m->n_nodes; i++) {
      free_ibnode(&m->nodes[i]);
    }
    free(m->nodes);
  }
  free(m);
}

/**
 * IBMeshManager
 */

typedef struct
{
  IBMesh *meshes;
  size_t n_meshes;
} IBMeshManager;

// #define IBMESH(i) (ibmeshmgr->meshes[i])

IBMeshManager *ibmeshmgr = NULL;

void init_ibmeshmanager(size_t n_meshes) {
  if (ibmeshmgr != NULL) return;
  ibmeshmgr = (IBMeshManager*) malloc(sizeof(IBMeshManager));
  ibmeshmgr->n_meshes = n_meshes;
  ibmeshmgr->meshes = (IBMesh*) calloc(n_meshes, sizeof(IBMesh));
  // for (size_t i = 0; i < n_meshes; i++) {
  //   init_ibmesh(&ibmeshmgr->meshes[i], 0);
  // }
}

void free_ibmeshmanager() {
  if (ibmeshmgr == NULL) return;
  for (size_t i = 0; i < ibmeshmgr->n_meshes; i++) {
    free_ibmesh(&ibmeshmgr->meshes[i]);
  }
  free(ibmeshmgr);
  ibmeshmgr = NULL;
}