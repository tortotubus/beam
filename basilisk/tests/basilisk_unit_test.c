
#include "tests/test.h"
#include "library/ibm/IBMempool.h"

// test_ibmempool.c
// #include <stddef.h>
// #include <stdint.h>
// #include <signal.h>

// #include <fenv.h>

// ----------------- helpers -----------------

static int ptr_in_active_(IBMempool *mp, IBNode *p)
{
  for (long i = 0; i < mp->na; i++)
    if (mp->active[i] == p) return 1;
  return 0;
}

static void (*g_death_fn_)(void) = NULL;

// ----------------- tests -----------------

static void test_init_sets_invariants(void)
{
  IBMempool mp = ibmempool_init(4096);

  TT_ASSERT_NONNULL(mp.pool);
  TT_ASSERT_PTR_EQ(mp.active, NULL);
  TT_ASSERT_EQ_LONG(mp.na, 0);
  TT_ASSERT_EQ_LONG(mp.nc, 0);

  ibmempool_free(&mp);

  TT_ASSERT_PTR_EQ(mp.pool, NULL);
  TT_ASSERT_PTR_EQ(mp.active, NULL);
  TT_ASSERT_EQ_LONG(mp.na, 0);
  TT_ASSERT_EQ_LONG(mp.nc, 0);
}

static void test_alloc_first_node(void)
{
  IBMempool mp = ibmempool_init(4096);

  IBNode *n0 = ibmempool_alloc_node(&mp);
  TT_ASSERT_NONNULL(n0);

  TT_ASSERT_NONNULL(mp.active);
  TT_ASSERT_EQ_LONG(mp.na, 1);
  TT_ASSERT_PTR_EQ(mp.active[0], n0);

  ibmempool_free(&mp);
}

static void test_alloc_many_unique_pointers(void)
{
  IBMempool mp = ibmempool_init(1u << 20);

  enum { N = 256 };
  IBNode *nodes[N];

  for (int i = 0; i < N; i++) {
    nodes[i] = ibmempool_alloc_node(&mp);
    TT_ASSERT_NONNULL(nodes[i]);
    TT_ASSERT_EQ_LONG(mp.na, i + 1);
  }

  // Sanity: check no duplicates among the first N allocations.
  // (If your mempool guarantees unique active allocations, this should hold.)
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      TT_ASSERT_PTR_NE(nodes[i], nodes[j]);

  ibmempool_free(&mp);
}

static void test_free_node_swap_with_last(void)
{
  IBMempool mp = ibmempool_init(1u << 20);

  IBNode *a = ibmempool_alloc_node(&mp);
  IBNode *b = ibmempool_alloc_node(&mp);
  IBNode *c = ibmempool_alloc_node(&mp);
  TT_ASSERT_EQ_LONG(mp.na, 3);

  // Remove middle (index 1). Should swap in last (c) and decrement na.
  ibmempool_free_node(&mp, 1);

  TT_ASSERT_EQ_LONG(mp.na, 2);
  TT_ASSERT_TRUE(ptr_in_active_(&mp, a));
  TT_ASSERT_FALSE(ptr_in_active_(&mp, b));
  TT_ASSERT_TRUE(ptr_in_active_(&mp, c));

  ibmempool_free(&mp);
}

static void test_free_node_last_element(void)
{
  IBMempool mp = ibmempool_init(1u << 20);

  IBNode *a = ibmempool_alloc_node(&mp);
  IBNode *b = ibmempool_alloc_node(&mp);
  TT_ASSERT_EQ_LONG(mp.na, 2);

  // Free last (index na-1)
  ibmempool_free_node(&mp, mp.na - 1);

  TT_ASSERT_EQ_LONG(mp.na, 1);
  TT_ASSERT_TRUE(ptr_in_active_(&mp, a));
  TT_ASSERT_FALSE(ptr_in_active_(&mp, b));

  ibmempool_free(&mp);
}

static IBMempool *g_mp_;
static long g_bad_index_;

static void call_free_node_bad_index_(void)
{
  ibmempool_free_node(g_mp_, g_bad_index_);
}


int main(void)
{
  TT_RUN(test_init_sets_invariants);
  TT_RUN(test_alloc_first_node);
  TT_RUN(test_alloc_many_unique_pointers);
  TT_RUN(test_free_node_swap_with_last);
  TT_RUN(test_free_node_last_element);

  return tt_summary_();
}
