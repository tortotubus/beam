#pragma once

#include "VTKHDF.hpp"
#include "VTKPolyData.hpp"

namespace beam {
class VTKHDFPolyData : public VTKHDF
{

protected:
  VTKPolyData& pd;

  static H5::DataSet create_and_write_connectivity_dataset(
    H5::Group& group,
    const std::vector<int64_t>& connectivity)
  {
    const hsize_t n = static_cast<hsize_t>(connectivity.size());

    // Filespace: rank-1, size = n, max = UNLIMITED
    hsize_t cur[1] = { n };
    hsize_t max[1] = { n };
    H5::DataSpace dspace(1, cur, max);

    // Must be chunked if UNLIMITED:
    // H5::DSetCreatPropList dcpl;
    // hsize_t chunk[1] = { std::max<hsize_t>( 1, std::min<hsize_t>(n ? n : 1024, 1024)) };
    // dcpl.setChunk(1, chunk);

    // Optional: better I/O / compression
    // dcpl.setShuffle();
    // dcpl.setDeflate(4);

    // File type is int64 little-endian
    H5::DataSet ds = group.createDataSet("Connectivity", H5::PredType::STD_I64LE, dspace);

    // Write entire dataset (use NATIVE type for portability)
    if (n > 0) { ds.write(connectivity.data(), H5::PredType::NATIVE_LLONG); }

    return ds;
  }

  static H5::DataSet create_and_write_offsets_dataset(
    H5::Group& group,
    const std::vector<int64_t>& offsets)
  {
    const hsize_t n = static_cast<hsize_t>(offsets.size());

    // Filespace: rank-1, size = n, max = UNLIMITED
    hsize_t cur[1] = { n };
    hsize_t max[1] = { n };
    H5::DataSpace dspace(1, cur, max);

    // Must be chunked if UNLIMITED:
    // H5::DSetCreatPropList dcpl;
    // hsize_t chunk[1] = { std::max<hsize_t>( 1, std::min<hsize_t>(n ? n : 1024, 1024)) };
    // dcpl.setChunk(1, chunk);

    // Optional: better I/O / compression
    // dcpl.setShuffle();
    // dcpl.setDeflate(4);

    // File type is int64 little-endian
    H5::DataSet ds = group.createDataSet("Offsets", H5::PredType::STD_I64LE, dspace);

    // Write entire dataset (use NATIVE type for portability)
    if (n > 0) { ds.write(offsets.data(), H5::PredType::NATIVE_LLONG); }

    return ds;
  }

  static H5::DataSet create_and_write_number_of_cells_dataset(
    H5::Group& group,
    const int64_t number_of_cells
  ) {
    hsize_t cur[1] = { 1 };
    hsize_t max[1] = { 1 };
    H5::DataSpace dspace(1, cur, max);
    H5::DataSet ds = group.createDataSet("NumberOfCells", H5::PredType::STD_I64LE, dspace);
    int64_t number_of_cells_data[1] = { number_of_cells };
    ds.write(number_of_cells_data, H5::PredType::NATIVE_LLONG);
    return ds;
  }

  static H5::DataSet create_and_write_number_of_connectivity_ids_dataset(
    H5::Group& group,
    const int64_t number_of_connectivity_ids
  ) {
    hsize_t cur[1] = { 1 };
    hsize_t max[1] = { 1 };
    H5::DataSpace dspace(1, cur, max);
    H5::DataSet ds = group.createDataSet("NumberOfConnectivityIds", H5::PredType::STD_I64LE, dspace);
    int64_t number_of_connectivity_ids_data[1] = { number_of_connectivity_ids };
    ds.write(number_of_connectivity_ids_data, H5::PredType::NATIVE_LLONG);
    return ds;
  }

  // H5::DataSet number_of_points;
  // H5::DataSet points;

  void write_points_datasets() {
    {
      const hsize_t n = static_cast<hsize_t>(pd.number_of_points());
      hsize_t cur[2] = { n , 3 };
      hsize_t max[2] = { n , 3 };
      H5::DataSpace dspace(2, cur, max);
      H5::DataSet ds = grp_vtkhdf.createDataSet("Points", H5::PredType::IEEE_F32LE, dspace);
      if (n > 0) { ds.write(pd.points.data(), H5::PredType::NATIVE_FLOAT); }
    }
    {
      int64_t data[1] = { static_cast<int64_t>(pd.number_of_points()) };
      hsize_t cur[1] = { 1 };
      hsize_t max[1] = { 1 };
      H5::DataSpace dspace(1, cur, max);
      H5::DataSet ds = grp_vtkhdf.createDataSet("NumberOfPoints", H5::PredType::STD_I64LE, dspace);
      ds.write(data, H5::PredType::NATIVE_LLONG);
    }
  }

  H5::Group vertices;

  // H5::DataSet vertices_connectivity;
  // H5::DataSet vertices_number_of_cells;
  // H5::DataSet vertices_number_of_connectivity_ids;
  // H5::DataSet vertices_offsets;

  void write_vertices_datasets() {     
    // Group: /VTKHDF/Lines
    vertices = ensure_group(file, "/VTKHDF/Vertices");
    // Dataset: /VTKHDF/Vertices/Connectivity
    create_and_write_connectivity_dataset(vertices, pd.vertices_connectivity);
    // Dataset: /VTKHDF/Vertices/NumberOfConnectivityIds
    create_and_write_number_of_connectivity_ids_dataset(vertices, pd.vertices_connectivity.size());
    // Dataset: /VTKHDF/Vertices/Offsets
    create_and_write_offsets_dataset(vertices, pd.vertices_offsets);
    // Dataset: /VTKHDF/Vertices/NumberOfCells
    create_and_write_number_of_cells_dataset(vertices, pd.number_of_vertices());
  }

  H5::Group lines;

  // H5::DataSet lines_connectivity;
  // H5::DataSet lines_number_of_cells;
  // H5::DataSet lines_number_of_connectivity_ids;
  // H5::DataSet lines_offsets;

  void write_lines_datasets() {     
    // Group: /VTKHDF/Lines
    lines = ensure_group(file, "/VTKHDF/Lines");
    // Dataset: /VTKHDF/Lines/Connectivity
    create_and_write_connectivity_dataset(lines, pd.lines_connectivity);
    // Dataset: /VTKHDF/Lines/NumberOfConnectivityIds
    create_and_write_number_of_connectivity_ids_dataset(lines, pd.lines_connectivity.size());
    // Dataset: /VTKHDF/Lines/Offsets
    create_and_write_offsets_dataset(lines, pd.lines_offsets);
    // Dataset: /VTKHDF/Lines/NumberOfCells
    create_and_write_number_of_cells_dataset(lines, pd.number_of_lines());
  }

  H5::Group strips;

  // H5::DataSet strips_connectivity;
  // H5::DataSet strips_number_of_cells;
  // H5::DataSet strips_number_of_connectivity_ids;
  // H5::DataSet strips_offsets;

  void write_strips_datasets() {     
    // Group: /VTKHDF/Strips
    strips = ensure_group(file, "/VTKHDF/Strips");
    // Dataset: /VTKHDF/Strips/Connectivity
    create_and_write_connectivity_dataset(strips, pd.strips_connectivity);
    // Dataset: /VTKHDF/Strips/NumberOfConnectivityIds
    create_and_write_number_of_connectivity_ids_dataset(strips, pd.strips_connectivity.size());
    // Dataset: /VTKHDF/Strips/Offsets
    create_and_write_offsets_dataset(strips, pd.strips_offsets);
    // Dataset: /VTKHDF/Strips/NumberOfCells
    create_and_write_number_of_cells_dataset(strips, pd.number_of_strips());
  }

  H5::Group polygons;

  // H5::DataSet polygons_connectivity;
  // H5::DataSet polygons_number_of_cells;
  // H5::DataSet polygons_number_of_connectivity_ids;
  // H5::DataSet polygons_offsets;

  void write_polygons_datasets() {     
    // Group: /VTKHDF/Polygons
    polygons = ensure_group(file, "/VTKHDF/Polygons");
    // Dataset: /VTKHDF/Polygons/Connectivity
    create_and_write_connectivity_dataset(polygons, pd.polygons_connectivity);
    // Dataset: /VTKHDF/Polygons/NumberOfConnectivityIds
    create_and_write_number_of_connectivity_ids_dataset(polygons, pd.polygons_connectivity.size());
    // Dataset: /VTKHDF/Polygons/Offsets
    create_and_write_offsets_dataset(polygons, pd.polygons_offsets);
    // Dataset: /VTKHDF/Polygons/NumberOfCells
    create_and_write_number_of_cells_dataset(polygons, pd.number_of_polygons());
  }
public:
  VTKHDFPolyData(std::string filename, VTKPolyData& pd)
    : VTKHDF(filename, POLYDATA)
    , pd(pd) {}; 

  void write()
  {
    write_points_datasets();
    write_vertices_datasets();
    write_lines_datasets();
    write_strips_datasets();
    write_polygons_datasets();
  }
};
}