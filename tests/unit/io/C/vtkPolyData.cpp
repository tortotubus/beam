// // vtkPolyDataTest.cpp
// #include <gtest/gtest.h>
// #include <cstdint>

// // Adjust include path / name as needed:
// #include "io/C/vtkPolyData.hpp"

// namespace ELFF {

// // Helper to get a clean, zero-initialized polydata
// static vtkPolyData make_empty_polydata()
// {
//   vtkPolyData pd{};
//   // enum values will be zero-initialized to BUILDING, but we can be explicit:
//   pd.points_state       = BUILDING;
//   pd.connectivity_state = BUILDING;
//   pd.fields_state       = BUILDING;
//   return pd;
// }

// // ------------------------- Points -------------------------

// TEST(VtkPolyDataTest, MallocAndFreePoints)
// {
//   vtkPolyData pd = make_empty_polydata();

//   EXPECT_EQ(pd.points, nullptr);
//   EXPECT_EQ(pd.n_points, 0u);
//   EXPECT_EQ(pd.m_points, 0u);

//   vtk_polydata_malloc_points(pd, 3); // capacity for 9 floats (3 points)

//   EXPECT_NE(pd.points, nullptr);
//   EXPECT_EQ(pd.n_points, 0u);
//   EXPECT_EQ(pd.m_points, 3u);

//   vtk_polydata_free_points(pd);

//   EXPECT_EQ(pd.points, nullptr);
//   EXPECT_EQ(pd.n_points, 0u);
//   EXPECT_EQ(pd.m_points, 0u);
// }

// TEST(VtkPolyDataTest, AddPointStoresCoordinatesAndIncrementsCount)
// {
//   vtkPolyData pd = make_empty_polydata();
//   vtk_polydata_malloc_points(pd, 9); // enough for 3 points

//   int64_t idx0 = vtk_polydata_add_point(pd, 1.f, 2.f, 3.f);
//   EXPECT_EQ(idx0, 0);
//   EXPECT_EQ(pd.n_points, 1u);
//   EXPECT_FLOAT_EQ(pd.points[0], 1.f);
//   EXPECT_FLOAT_EQ(pd.points[1], 2.f);
//   EXPECT_FLOAT_EQ(pd.points[2], 3.f);

//   // These expectations describe the *intended* semantics.
//   // With your current implementation, vtk_polydata_number_of_points()
//   // returns pd.n_points / 3, so this test will fail and highlight that bug.
//   EXPECT_EQ(vtk_polydata_number_of_points(pd), 1u);

//   vtk_polydata_free_points(pd);
// }

// // ---------------------- Vertices --------------------------

// TEST(VtkPolyDataTest, MallocAndFreeVertices)
// {
//   vtkPolyData pd = make_empty_polydata();

//   EXPECT_EQ(pd.vertices_connectivity, nullptr);
//   EXPECT_EQ(pd.vertices_offsets, nullptr);

//   vtk_polydata_malloc_vertices(pd, 4);

//   EXPECT_NE(pd.vertices_connectivity, nullptr);
//   EXPECT_NE(pd.vertices_offsets, nullptr);

//   EXPECT_EQ(pd.m_vertices_connectivity, 4u);
//   EXPECT_EQ(pd.n_vertices_connectivity, 0u);

//   EXPECT_EQ(pd.m_vertices_offsets, 5u);  // n + 1
//   EXPECT_EQ(pd.n_vertices_offsets, 1u);  // only offset[0] valid
//   EXPECT_EQ(pd.vertices_offsets[0], 0);

//   vtk_polydata_free_vertices(pd);

//   EXPECT_EQ(pd.vertices_connectivity, nullptr);
//   EXPECT_EQ(pd.vertices_offsets, nullptr);
//   EXPECT_EQ(pd.m_vertices_connectivity, 0u);
//   EXPECT_EQ(pd.n_vertices_connectivity, 0u);
//   EXPECT_EQ(pd.m_vertices_offsets, 0u);
//   EXPECT_EQ(pd.n_vertices_offsets, 0u);
// }

// TEST(VtkPolyDataTest, AddVertexUsesExistingPointAndUpdatesOffsets)
// {
//   vtkPolyData pd = make_empty_polydata();

//   // Set up one point
//   vtk_polydata_malloc_points(pd, 3); // capacity for 1 point (3 floats)
//   vtk_polydata_add_point(pd, 0.f, 0.f, 0.f);

//   // Allocate room for one vertex
//   vtk_polydata_malloc_vertices(pd, 1);

//   int64_t v_id = vtk_polydata_add_vertex(pd, 0);

//   EXPECT_EQ(v_id, 0);
//   EXPECT_EQ(pd.n_vertices_connectivity, 1u);
//   EXPECT_EQ(pd.vertices_connectivity[0], 0);

//   EXPECT_EQ(pd.n_vertices_offsets, 2u);
//   EXPECT_EQ(pd.vertices_offsets[0], 0);
//   EXPECT_EQ(pd.vertices_offsets[1], 1);

//   EXPECT_EQ(vtk_polydata_number_of_vertices(pd), 1u);
// }

// // This requires gtest's death test support and that you link with gtest_main.
// TEST(VtkPolyDataDeathTest, AddVertexWithInvalidPointAborts)
// {
//   vtkPolyData pd = make_empty_polydata();

//   // No points, but allocate vertices
//   vtk_polydata_malloc_vertices(pd, 1);

//   EXPECT_DEATH(
//       {
//         vtk_polydata_add_vertex(pd, 0);  // point 0 does not exist
//       },
//       ".*");
// }

// // ------------------------ Lines ---------------------------

// TEST(VtkPolyDataTest, MallocAndFreeLines)
// {
//   vtkPolyData pd = make_empty_polydata();

//   vtk_polydata_malloc_lines(pd, 4);

//   EXPECT_NE(pd.lines_connectivity, nullptr);
//   EXPECT_NE(pd.lines_offsets, nullptr);

//   EXPECT_EQ(pd.m_lines_connectivity, 4u);
//   EXPECT_EQ(pd.n_lines_connectivity, 0u);

//   EXPECT_EQ(pd.m_lines_offsets, 5u);
//   EXPECT_EQ(pd.n_lines_offsets, 1u);
//   EXPECT_EQ(pd.lines_offsets[0], 0);

//   vtk_polydata_free_lines(pd);

//   EXPECT_EQ(pd.lines_connectivity, nullptr);
//   EXPECT_EQ(pd.lines_offsets, nullptr);
//   EXPECT_EQ(pd.m_lines_connectivity, 0u);
//   EXPECT_EQ(pd.n_lines_connectivity, 0u);
//   EXPECT_EQ(pd.m_lines_offsets, 0u);
//   EXPECT_EQ(pd.n_lines_offsets, 0u);
// }

// TEST(VtkPolyDataTest, AddLineBetweenTwoValidPoints)
// {
//   vtkPolyData pd = make_empty_polydata();

//   // Two points
//   vtk_polydata_malloc_points(pd, 6); // capacity for 2 points / 6 floats
//   vtk_polydata_add_point(pd, 0.f, 0.f, 0.f); // index 0
//   vtk_polydata_add_point(pd, 1.f, 0.f, 0.f); // index 1

//   // Allocate lines
//   vtk_polydata_malloc_lines(pd, 2); // room for a single 2-point line

//   int64_t line_id = vtk_polydata_add_line(pd, 0, 1);

//   EXPECT_EQ(line_id, 0);

//   EXPECT_EQ(pd.n_lines_connectivity, 2u);
//   EXPECT_EQ(pd.lines_connectivity[0], 0);
//   EXPECT_EQ(pd.lines_connectivity[1], 1);

//   EXPECT_EQ(pd.lines_offsets[0], 0);
//   EXPECT_EQ(pd.lines_offsets[1], 2);

//   // Again, this expectation encodes the *intended* behavior.
//   // Your current implementation increments n_vertices_offsets instead
//   // of n_lines_offsets, so n_lines_offsets stays 1 and this will fail,
//   // flagging that bug.
//   EXPECT_EQ(pd.n_lines_offsets, 2u);
//   EXPECT_EQ(vtk_polydata_number_of_lines(pd), 1u);
// }

// TEST(VtkPolyDataDeathTest, AddLineWithInvalidPointAborts)
// {
//   vtkPolyData pd = make_empty_polydata();

//   vtk_polydata_malloc_points(pd, 3);
//   vtk_polydata_add_point(pd, 0.f, 0.f, 0.f); // only point 0 exists

//   vtk_polydata_malloc_lines(pd, 2);

//   EXPECT_DEATH(
//       {
//         vtk_polydata_add_line(pd, 0, 1); // point 1 does not exist
//       },
//       ".*");
// }

// // ---------------------- Sealing logic ---------------------

// TEST(VtkPolyDataDeathTest, AddingPointsAfterConnectivityAborts)
// {
//   vtkPolyData pd = make_empty_polydata();

//   vtk_polydata_malloc_points(pd, 6);
//   vtk_polydata_add_point(pd, 0.f, 0.f, 0.f);

//   vtk_polydata_malloc_vertices(pd, 1);
//   vtk_polydata_add_vertex(pd, 0); // this seals points_state

//   EXPECT_TRUE(vtk_polydata_points_is_sealed(pd));
//   EXPECT_FALSE(vtk_polydata_connectivity_is_sealed(pd));

//   // After connectivity has been added, further point additions should abort
//   EXPECT_DEATH(
//       {
//         vtk_polydata_add_point(pd, 1.f, 1.f, 1.f);
//       },
//       ".*");
// }

// } // namespace ELFF
