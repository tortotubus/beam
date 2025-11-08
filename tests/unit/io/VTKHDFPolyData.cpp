#include <gtest/gtest.h>
#include <H5Cpp.h>
#include <filesystem>
#include <string>

// Your headers
#include "io/VTKPolyData.hpp"
#include "io/VTKHDFPolyData.hpp"

using namespace H5;
using namespace beam;

namespace {

// ---- tiny HDF5 helpers ----------------------------------------------------

bool linkExists(H5::H5File& f, const std::string& path) {
  return H5Lexists(f.getId(), path.c_str(), H5P_DEFAULT) > 0;
}

hsize_t oneDDatasetLength(const H5::DataSet& ds) {
  H5::DataSpace sp = ds.getSpace();
  int rank = sp.getSimpleExtentNdims();
  EXPECT_EQ(rank, 1);
  hsize_t dims[1] = {0};
  sp.getSimpleExtentDims(dims, nullptr);
  return dims[0];
}

long long readScalarI64(const H5::DataSet& ds) {
  long long v = -1;
  H5::DataSpace fsp = ds.getSpace();
  hsize_t one[1] = {1};
  H5::DataSpace msp(1, one);
  ds.read(&v, H5::PredType::NATIVE_LLONG, msp, fsp);
  return v;
}

void expectSectionEmpty(H5::H5File& f, const char* sectionPath) {
  // e.g. "/VTKHDF/Vertices"
  ASSERT_TRUE(linkExists(f, sectionPath)) << sectionPath << " missing";

  std::string base(sectionPath);
  auto open = [&](const char* name) {
    return f.openDataSet(base + "/" + name);
  };

  // Connectivity: length 0
  {
    ASSERT_TRUE(linkExists(f, base + "/Connectivity"));
    auto ds = open("Connectivity");
    EXPECT_EQ(oneDDatasetLength(ds), 0u);
  }
  // Offsets: at least one entry (should be [0] when empty)
  {
    ASSERT_TRUE(linkExists(f, base + "/Offsets"));
    auto ds = open("Offsets");
    hsize_t n = oneDDatasetLength(ds);
    ASSERT_GE(n, 1u);
    // Read first value
    std::vector<long long> buf(n, -1);
    H5::DataSpace fsp = ds.getSpace();
    H5::DataSpace msp(1, &n);
    ds.read(buf.data(), H5::PredType::NATIVE_LLONG, msp, fsp);
    EXPECT_EQ(buf[0], 0) << sectionPath << "/Offsets[0] must be 0";
  }
  // NumberOfCells == 0
  {
    ASSERT_TRUE(linkExists(f, base + "/NumberOfCells"));
    auto ds = open("NumberOfCells");
    EXPECT_EQ(oneDDatasetLength(ds), 1u);
    EXPECT_EQ(readScalarI64(ds), 0);
  }
  // NumberOfConnectivityIds == 0
  {
    ASSERT_TRUE(linkExists(f, base + "/NumberOfConnectivityIds"));
    auto ds = open("NumberOfConnectivityIds");
    EXPECT_EQ(oneDDatasetLength(ds), 1u);
    EXPECT_EQ(readScalarI64(ds), 0);
  }
}

} // namespace

// ---- the test --------------------------------------------------------------

TEST(VTKHDFPolyData, WriteEmptyPolyDataDoesNotThrowAndLooksSane) {
  // 1) Build a default/empty polydata
  VTKPolyData pd;

  // If your VTKPolyData default ctor doesn't seed offsets with {0}, ensure it.
  // if (pd.vertices_offsets.empty())  pd.vertices_offsets.push_back(0);
  // if (pd.lines_offsets.empty())     pd.lines_offsets.push_back(0);
  // if (pd.strips_offsets.empty())    pd.strips_offsets.push_back(0);
  // if (pd.polygons_offsets.empty())  pd.polygons_offsets.push_back(0);

  // 2) Write it
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "empty_polydata.hdf";
  {
    VTKHDFPolyData writer(tmp.string(), pd);
    ASSERT_NO_THROW(writer.write());
  }

  // 3) Re-open and sanity-check structure
  H5::H5File file(tmp.string(), H5F_ACC_RDONLY);

  // Header group exists
  ASSERT_TRUE(linkExists(file, "/VTKHDF"));

  // Type attribute exists (optional: check equals "PolyData")
  {
    auto g = file.openGroup("/VTKHDF");
    ASSERT_TRUE(H5Aexists(g.getId(), "Type") > 0);
    // Read into a small buffer without assuming exact fixed length
    H5::Attribute a = g.openAttribute("Type");
    H5::DataType  t = a.getDataType();
    // Make a buffer of up to, say, 32 chars
    std::vector<char> buf(32, 0);
    // If it's fixed-length, read as-is; if it were variable, you'd handle differently.
    a.read(t, buf.data());
    // We won't ASSERT the exact string length; a basic sanity check:
    // EXPECT_STREQ("PolyData", buf.data());
  }

  // NumberOfPoints == 0 and Points dataset exists with (0,3)
  {
    ASSERT_TRUE(linkExists(file, "/VTKHDF/NumberOfPoints"));
    auto dsN = file.openDataSet("/VTKHDF/NumberOfPoints");
    EXPECT_EQ(oneDDatasetLength(dsN), 1u);
    EXPECT_EQ(readScalarI64(dsN), 0);

    ASSERT_TRUE(linkExists(file, "/VTKHDF/Points"));
    auto dsP = file.openDataSet("/VTKHDF/Points");
    H5::DataSpace sp = dsP.getSpace();
    int rank = sp.getSimpleExtentNdims();
    EXPECT_EQ(rank, 2);
    hsize_t dims[2] = {0,0};
    sp.getSimpleExtentDims(dims, nullptr);
    EXPECT_EQ(dims[0], 0u); // zero rows
    EXPECT_EQ(dims[1], 3u); // xyz
  }

  // Each cell section exists and is empty but well-formed
  expectSectionEmpty(file, "/VTKHDF/Vertices");
  expectSectionEmpty(file, "/VTKHDF/Lines");
  expectSectionEmpty(file, "/VTKHDF/Strips");
  expectSectionEmpty(file, "/VTKHDF/Polygons");
}
