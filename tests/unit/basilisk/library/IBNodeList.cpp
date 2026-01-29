// test_ibnodelist.cc
#include <gtest/gtest.h>
#include <cstddef>
#include <cstring>

extern "C" {

    typedef struct {
      double x, y, z;
    } coord;

    typedef int Index;

    typedef struct {
      size_t n;    /*!< Current number of entries in the stencil */
      size_t nm;   /*!< Allocated capacity of the stencil */
      Index  *p;   /*!< Pointer to stencil entries (owned) */
    } Cache;

    #define PESKIN_SUPPORT_RADIUS 2
    #define sq(x) ((x)*(x))
    #define cube(x) ((x)*(x)*(x))

    #include "library/ibm/IBNode.h"
    #include "library/ibm/IBNodeList.h"
}

namespace {

static IBNode MakeNode(double base) {
  IBNode n;
  ibnode_init(&n);

  // Fill a few fields with recognizable values.
  // (coord typically has x,y,z even if dimension<3.)
  n.lagpos.x = base + 1.0;
  n.lagpos.y = base + 2.0;
  n.lagpos.z = base + 3.0;

  n.force.x  = base + 4.0;
  n.force.y  = base + 5.0;
  n.force.z  = base + 6.0;

  n.weight   = base + 7.0;
  n.pid      = static_cast<int>(base);

  // If your IBNode has an owned stencil buffer like earlier, seed it too.
  // This makes deep-copy tests meaningful.
  if (n.stencil.p && n.stencil.nm > 0) {
    n.stencil.n = 1;
    n.stencil.p[0] = n.stencil.p[0]; // no-op; keep compiler quiet if Index isn't POD
  }

  return n;
}

static void ExpectNodeEq(const IBNode& a, const IBNode& b) {
  EXPECT_DOUBLE_EQ(a.lagpos.x, b.lagpos.x);
  EXPECT_DOUBLE_EQ(a.lagpos.y, b.lagpos.y);
  EXPECT_DOUBLE_EQ(a.lagpos.z, b.lagpos.z);

  EXPECT_DOUBLE_EQ(a.force.x,  b.force.x);
  EXPECT_DOUBLE_EQ(a.force.y,  b.force.y);
  EXPECT_DOUBLE_EQ(a.force.z,  b.force.z);

  EXPECT_DOUBLE_EQ(a.weight, b.weight);
  EXPECT_EQ(a.pid, b.pid);
}

} // namespace

TEST(IBNodeList, InitNullPointerFails) {
  EXPECT_EQ(ibnodelist_init(nullptr, 4), -1);
}

TEST(IBNodeList, InitZeroCapacityIsValidEmptyList) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 0), 0);
  EXPECT_EQ(L.size, 0u);
  EXPECT_EQ(L.capacity, 0u);
  EXPECT_EQ(L.nodes, nullptr);

  // Should be safe / no-op.
  ibnodelist_free(&L);
  EXPECT_EQ(L.nodes, nullptr);
  EXPECT_EQ(L.size, 0u);
  EXPECT_EQ(L.capacity, 0u);
}

TEST(IBNodeList, ReserveIncreasesCapacity) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 0), 0);

  ASSERT_EQ(ibnodelist_reserve(&L, 4), 0);
  EXPECT_GE(L.capacity, 4u);
  EXPECT_NE(L.nodes, nullptr);

  ibnodelist_free(&L);
}

TEST(IBNodeList, GetOutOfRangeReturnsNull) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 0), 0);

  EXPECT_EQ(ibnodelist_get(&L, 0), nullptr);
  EXPECT_EQ(ibnodelist_get_const(&L, 0), nullptr);

  IBNode n = MakeNode(10.0);
  ASSERT_GE(ibnodelist_add_copy(&L, &n), 0);

  EXPECT_EQ(ibnodelist_get(&L, 1), nullptr);
  EXPECT_EQ(ibnodelist_get_const(&L, 1), nullptr);

  ibnode_free(&n);
  ibnodelist_free(&L);
}

TEST(IBNodeList, AddCopyFromZeroCapacityGrowsAndStoresValue) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 0), 0);

  IBNode src = MakeNode(1.0);
  int idx = ibnodelist_add_copy(&L, &src);
  ASSERT_EQ(idx, 0);
  ASSERT_EQ(L.size, 1u);
  ASSERT_GE(L.capacity, 1u);

  const IBNode* got = ibnodelist_get_const(&L, 0);
  ASSERT_NE(got, nullptr);
  ExpectNodeEq(*got, src);

  ibnode_free(&src);
  ibnodelist_free(&L);
}

TEST(IBNodeList, AddCopyIsDeepCopy_SourceMutationDoesNotAffectStoredNode) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 0), 0);

  IBNode src = MakeNode(2.0);
  ASSERT_GE(ibnodelist_add_copy(&L, &src), 0);

  IBNode* stored = ibnodelist_get(&L, 0);
  ASSERT_NE(stored, nullptr);

  // Keep original values for comparison
  IBNode before = *stored;

  // Mutate source after insertion
  src.lagpos.x += 1000.0;
  src.force.y  -= 1000.0;
  src.weight   += 1000.0;
  src.pid      += 1000;

  // If stencil is owned/dynamic, mutate source buffer too.
  // Deep copy should ensure stored->stencil.p is a different allocation (or at least independent).
  if (src.stencil.p && stored->stencil.p && src.stencil.nm > 0) {
    // This assertion will FAIL if ibnode_copy is shallow (good: it catches the bug).
    EXPECT_NE(src.stencil.p, stored->stencil.p);

    // Try to mutate source buffer content (only safe if Index is trivially assignable).
    // If Index isn't trivially assignable, remove this part and keep pointer inequality check.
    src.stencil.n = 0; // benign mutation
  }

  // Stored node should remain unchanged compared to snapshot taken before mutating src
  ExpectNodeEq(*stored, before);

  ibnode_free(&src);
  ibnodelist_free(&L);
}

TEST(IBNodeList, FreeResetsListState) {
  IBNodeList L{};
  ASSERT_EQ(ibnodelist_init(&L, 2), 0);

  IBNode src = MakeNode(3.0);
  ASSERT_GE(ibnodelist_add_copy(&L, &src), 0);
  ibnode_free(&src);

  ibnodelist_free(&L);

  EXPECT_EQ(L.nodes, nullptr);
  EXPECT_EQ(L.size, 0u);
  EXPECT_EQ(L.capacity, 0u);
}
