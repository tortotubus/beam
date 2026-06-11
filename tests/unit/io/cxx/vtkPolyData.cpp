#include <gtest/gtest.h>

#include <elff/io/C/vtkPolyData.hpp>
#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/io/CXX/vtkPolyData.hpp>

#include <hdf5.h>

#include <filesystem>
#include <string>
#include <vector>

namespace ELFF {
namespace IO {
namespace CXX {
namespace {

class InspectablePolyData : public vtkPolyData {
public:
  InspectablePolyData() = default;
  explicit InspectablePolyData(const vtkPolyData& other)
    : vtkPolyData(other)
  {
  }

  using vtkPolyData::to_c_struct;
};

std::filesystem::path
test_output_path(const char* name)
{
  return std::filesystem::path(name);
}

void
remove_file_if_present(const std::filesystem::path& path)
{
  std::error_code ec;
  std::filesystem::remove(path, ec);
}

vtkPolyData
make_triangle_polydata(double value_offset)
{
  vtkPolyData pd;

  const int64_t p0 = pd.add_point(static_cast<float>(value_offset), 0.f, 0.f);
  const int64_t p1 = pd.add_point(static_cast<float>(value_offset + 1.), 0.f, 0.f);
  const int64_t p2 = pd.add_point(static_cast<float>(value_offset), 1.f, 0.f);
  pd.add_polygon({ p0, p1, p2 });

  const int64_t point_id = pd.add_pointdata_scalar("point_id");
  std::vector<double>& point_data = pd.get_pointdata(point_id);
  point_data[0] = value_offset + 0.;
  point_data[1] = value_offset + 1.;
  point_data[2] = value_offset + 2.;

  const int64_t area_id = pd.add_celldata_scalar("triangle_area");
  std::vector<double>& area = pd.get_celldata(area_id);
  area[0] = value_offset + 0.5;

  return pd;
}

bool
hdf5_link_exists(const std::filesystem::path& path, const char* link_path)
{
  hid_t file_id = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return false;

  const htri_t exists = H5Lexists(file_id, link_path, H5P_DEFAULT);
  H5Fclose(file_id);
  return exists > 0;
}

TEST(vtkPolyDataCXXTest, AddsCellDataFields)
{
  vtkPolyData pd;

  const int64_t p0 = pd.add_point(0.f, 0.f, 0.f);
  const int64_t p1 = pd.add_point(1.f, 0.f, 0.f);
  const int64_t p2 = pd.add_point(1.f, 1.f, 0.f);
  const int64_t p3 = pd.add_point(0.f, 1.f, 0.f);

  pd.add_vertex(p0);
  pd.add_line(p1, p2);
  pd.add_polygon({ p0, p2, p3 });

  ASSERT_EQ(pd.number_of_cells(), 3);

  const int64_t scalar_id = pd.add_celldata_scalar("cell_id");
  const int64_t vector_id = pd.add_celldata_vector("cell_vec", 3);

  std::vector<double>& scalar = pd.get_celldata(scalar_id);
  scalar[0] = 10.;
  scalar[1] = 11.;
  scalar[2] = 12.;

  pd.set_celldata_vector3(vector_id, 0, { 1., 0., 0. });
  pd.set_celldata_vector3(vector_id, 1, { 0., 1., 0. });
  pd.set_celldata_vector3(vector_id, 2, { 0., 0., 1. });

  ASSERT_EQ(pd.number_of_celldata(), 2);
  EXPECT_EQ(pd.get_celldata(scalar_id).size(), 3);
  EXPECT_EQ(pd.get_celldata(vector_id).size(), 9);
  EXPECT_DOUBLE_EQ(pd.get_celldata(scalar_id)[2], 12.);
  EXPECT_DOUBLE_EQ(pd.get_celldata(vector_id)[8], 1.);
}

TEST(vtkPolyDataCXXTest, ConvertsCellDataToCStruct)
{
  InspectablePolyData pd;

  const int64_t p0 = pd.add_point(0.f, 0.f, 0.f);
  const int64_t p1 = pd.add_point(1.f, 0.f, 0.f);
  const int64_t p2 = pd.add_point(0.f, 1.f, 0.f);
  pd.add_polygon({ p0, p1, p2 });

  const int64_t cell_id = pd.add_celldata_vector("triangle_normal", 3);
  pd.set_celldata_vector3(cell_id, 0, { 0., 0., 1. });

  C::vtkPolyData c_pd = pd.to_c_struct();

  ASSERT_EQ(c_pd.n_celldata, 1);
  ASSERT_STREQ(c_pd.celldata_names[0], "triangle_normal");
  ASSERT_EQ(c_pd.celldata_ncomp[0], 3);
  EXPECT_DOUBLE_EQ(c_pd.celldata_data[0][0], 0.);
  EXPECT_DOUBLE_EQ(c_pd.celldata_data[0][1], 0.);
  EXPECT_DOUBLE_EQ(c_pd.celldata_data[0][2], 1.);

  C::vtk_polydata_free(&c_pd);
}

TEST(vtkPolyDataCXXTest, AppendManyPreservesCellDataAndOffsetsConnectivity)
{
  vtkPolyData first = make_triangle_polydata(0.);
  vtkPolyData second = make_triangle_polydata(10.);

  vtkPolyData appended = vtkPolyData::append_many({ first, second });

  ASSERT_EQ(appended.number_of_points(), 6);
  ASSERT_EQ(appended.number_of_polygons(), 2);
  ASSERT_EQ(appended.number_of_pointdata(), 1);
  ASSERT_EQ(appended.number_of_celldata(), 1);

  const std::vector<double>& point_data = appended.get_pointdata(0);
  const std::vector<double>& cell_data = appended.get_celldata(0);

  ASSERT_EQ(point_data.size(), 6);
  ASSERT_EQ(cell_data.size(), 2);
  EXPECT_DOUBLE_EQ(point_data[0], 0.);
  EXPECT_DOUBLE_EQ(point_data[3], 10.);
  EXPECT_DOUBLE_EQ(cell_data[0], 0.5);
  EXPECT_DOUBLE_EQ(cell_data[1], 10.5);

  InspectablePolyData inspectable(appended);
  C::vtkPolyData c_pd = inspectable.to_c_struct();

  ASSERT_EQ(c_pd.n_polygons_offsets, 3);
  ASSERT_EQ(c_pd.n_polygons_connectivity, 6);
  EXPECT_EQ(c_pd.polygons_connectivity[3], 3);
  EXPECT_EQ(c_pd.polygons_connectivity[4], 4);
  EXPECT_EQ(c_pd.polygons_connectivity[5], 5);

  C::vtk_polydata_free(&c_pd);
}

TEST(vtkPolyDataCXXTest, WritesStaticHDFWithCellData)
{
  const std::filesystem::path output_path =
    test_output_path("elff_vtk_polydata_cxx_static_fields.vtkhdf");
  remove_file_if_present(output_path);

  vtkPolyData pd = make_triangle_polydata(0.);

  vtkHDFPolyData writer(output_path.string(), pd);
  writer.write_new_static(true);

  ASSERT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);
  EXPECT_TRUE(hdf5_link_exists(output_path, "/VTKHDF/PointData/point_id"));
  EXPECT_TRUE(hdf5_link_exists(output_path, "/VTKHDF/CellData/triangle_area"));
}

} // namespace
} // namespace CXX
} // namespace IO
} // namespace ELFF
