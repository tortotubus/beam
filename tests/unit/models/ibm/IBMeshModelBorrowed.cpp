#include <gtest/gtest.h>

#include "elff/c/models/beam/IBEulerBeamHuang.h"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/ibm/IBRuntime.h"

#include <cstdlib>

#include "basilisk/library/ibm/IBMeshModel.h"

namespace {

TEST(IBMeshModelBorrowedTest, DestroyOnlyFreesAdapterTables)
{
  ib_runtime_t runtime = ib_runtime_new();
  ASSERT_NE(runtime, nullptr);
  ib_euler_beam_bcs_t bcs = {};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_FREE;
  const vertex_t origin = {0, 0, 0};
  ib_euler_beam_huang_t beam = ib_euler_beam_huang_new(
    runtime, origin, bcs, 1.0, 0.01, 1.0, 6);
  ASSERT_NE(beam, nullptr);

  IBMeshModel adapter = ibmeshmodel_force_coupled_init();
  adapter.ctx = beam;
  adapter.force_ops->destroy = nullptr;
  ibmeshmodel_destroy(&adapter);

  EXPECT_EQ(adapter.ctx, nullptr);
  EXPECT_EQ(adapter.force_ops, nullptr);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(beam), 6);
  ib_runtime_delete(runtime);
}

} // namespace
