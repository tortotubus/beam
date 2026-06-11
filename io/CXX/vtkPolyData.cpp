#include "elff/io/CXX/vtkPolyData.hpp"

#include <cstdlib>   // malloc, free
#include <cstring>   // std::memcpy

namespace ELFF {
namespace IO {
namespace CXX {

C::vtkPolyData
vtkPolyData::to_c_struct() const
{
  // Points:
  const size_t points_size = this->points.size();
  const size_t n_points = this->number_of_points();

  float* points = static_cast<float*>(std::malloc(sizeof(float) * points_size));

  if (!points) {
    ELFF_ABORT("malloc failure for points\n");
  }

  if (points_size > 0) {
    std::memcpy(points, this->points.data(), points_size * sizeof(float));
  }

  // Vertices: connectivity
  const size_t n_vertices_connectivity = this->vertices_connectivity.size();

  int64_t* vertices_connectivity = static_cast<int64_t*>(
    std::malloc(sizeof(int64_t) * n_vertices_connectivity));

  if (!vertices_connectivity) {
    ELFF_ABORT("malloc failure for vertices_connectivity\n");
  }

  if (n_vertices_connectivity > 0) {
    std::memcpy(vertices_connectivity,
                this->vertices_connectivity.data(),
                n_vertices_connectivity * sizeof(int64_t));
  }

  // Vertices: offsets
  const size_t n_vertices_offsets = this->vertices_offsets.size();

  int64_t* vertices_offsets =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_vertices_offsets));

  if (!vertices_offsets) {
    ELFF_ABORT("malloc failure for vertices_offsets\n");
  }

  if (n_vertices_offsets > 0) {
    std::memcpy(vertices_offsets,
                this->vertices_offsets.data(),
                n_vertices_offsets * sizeof(int64_t));
  }

  // Lines: connectivity
  const size_t n_lines_connectivity = this->lines_connectivity.size();

  int64_t* lines_connectivity =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_lines_connectivity));

  if (!lines_connectivity) {
    ELFF_ABORT("malloc failure for lines_connectivity\n");
  }

  if (n_lines_connectivity > 0) {
    std::memcpy(lines_connectivity,
                this->lines_connectivity.data(),
                n_lines_connectivity * sizeof(int64_t));
  }

  // Lines: offsets
  const size_t n_lines_offsets = this->lines_offsets.size();

  int64_t* lines_offsets =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_lines_offsets));

  if (!lines_offsets) {
    ELFF_ABORT("malloc failure for lines_offsets\n");
  }

  if (n_lines_offsets > 0) {
    std::memcpy(lines_offsets,
                this->lines_offsets.data(),
                n_lines_offsets * sizeof(int64_t));
  }

  // Strips: connectivity
  const size_t n_strips_connectivity = this->strips_connectivity.size();

  int64_t* strips_connectivity =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_strips_connectivity));

  if (!strips_connectivity) {
    ELFF_ABORT("malloc failure for strips_connectivity\n");
  }

  if (n_strips_connectivity > 0) {
    std::memcpy(strips_connectivity,
                this->strips_connectivity.data(),
                n_strips_connectivity * sizeof(int64_t));
  }

  // Strips: offsets
  const size_t n_strips_offsets = this->strips_offsets.size();

  int64_t* strips_offsets =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_strips_offsets));

  if (!strips_offsets) {
    ELFF_ABORT("malloc failure for strips_offsets\n");
  }

  if (n_strips_offsets > 0) {
    std::memcpy(strips_offsets,
                this->strips_offsets.data(),
                n_strips_offsets * sizeof(int64_t));
  }

  // Polygons: connectivity
  const size_t n_polygons_connectivity = this->polygons_connectivity.size();

  int64_t* polygons_connectivity = static_cast<int64_t*>(
    std::malloc(sizeof(int64_t) * n_polygons_connectivity));

  if (!polygons_connectivity) {
    ELFF_ABORT("malloc failure for polygons_connectivity\n");
  }

  if (n_polygons_connectivity > 0) {
    std::memcpy(polygons_connectivity,
                this->polygons_connectivity.data(),
                n_polygons_connectivity * sizeof(int64_t));
  }

  // Polygons: offsets
  const size_t n_polygons_offsets = this->polygons_offsets.size();

  int64_t* polygons_offsets =
    static_cast<int64_t*>(std::malloc(sizeof(int64_t) * n_polygons_offsets));

  if (!polygons_offsets) {
    ELFF_ABORT("malloc failure for polygons_offsets\n");
  }

  if (n_polygons_offsets > 0) {
    std::memcpy(polygons_offsets,
                this->polygons_offsets.data(),
                n_polygons_offsets * sizeof(int64_t));
  }

  const size_t n_pointdata = pointdata_data.size();

  char** pointdata_names = static_cast<char**>(
    std::calloc(n_pointdata, sizeof(char*)));
  size_t* pointdata_ncomp = static_cast<size_t*>(
    std::calloc(n_pointdata, sizeof(size_t)));
  double** pointdata_data = static_cast<double**>(
    std::calloc(n_pointdata, sizeof(double*)));

  if ((n_pointdata > 0) &&
      ((pointdata_names == nullptr) || (pointdata_ncomp == nullptr) ||
       (pointdata_data == nullptr))) {
    ELFF_ABORT("malloc failure for pointdata\n");
  }

  for (size_t i = 0; i < n_pointdata; ++i) {
    const std::string& field_name = this->pointdata_names[i];
    pointdata_names[i] = static_cast<char*>(
      std::malloc((field_name.size() + 1) * sizeof(char)));

    if (!pointdata_names[i]) {
      ELFF_ABORT("malloc failure for pointdata_names[i]\n");
    }

    std::memcpy(pointdata_names[i], field_name.c_str(), field_name.size() + 1);
    pointdata_ncomp[i] = this->pointdata_ncomp[i];

    const size_t field_size = this->pointdata_data[i].size();
    pointdata_data[i] = static_cast<double*>(
      std::malloc(field_size * sizeof(double)));

    if ((field_size > 0) && !pointdata_data[i]) {
      ELFF_ABORT("malloc failure for pointdata_data[i]\n");
    }

    if (field_size > 0) {
      std::memcpy(pointdata_data[i],
                  this->pointdata_data[i].data(),
                  field_size * sizeof(double));
    }
  }

  const size_t n_celldata = celldata_data.size();

  char** celldata_names = static_cast<char**>(
    std::calloc(n_celldata, sizeof(char*)));
  size_t* celldata_ncomp = static_cast<size_t*>(
    std::calloc(n_celldata, sizeof(size_t)));
  double** celldata_data = static_cast<double**>(
    std::calloc(n_celldata, sizeof(double*)));

  if ((n_celldata > 0) &&
      ((celldata_names == nullptr) || (celldata_ncomp == nullptr) ||
       (celldata_data == nullptr))) {
    ELFF_ABORT("malloc failure for celldata\n");
  }

  for (size_t i = 0; i < n_celldata; ++i) {
    const std::string& field_name = this->celldata_names[i];
    celldata_names[i] = static_cast<char*>(
      std::malloc((field_name.size() + 1) * sizeof(char)));

    if (!celldata_names[i]) {
      ELFF_ABORT("malloc failure for celldata_names[i]\n");
    }

    std::memcpy(celldata_names[i], field_name.c_str(), field_name.size() + 1);
    celldata_ncomp[i] = this->celldata_ncomp[i];

    const size_t field_size = this->celldata_data[i].size();
    celldata_data[i] = static_cast<double*>(
      std::malloc(field_size * sizeof(double)));

    if ((field_size > 0) && !celldata_data[i]) {
      ELFF_ABORT("malloc failure for celldata_data[i]\n");
    }

    if (field_size > 0) {
      std::memcpy(celldata_data[i],
                  this->celldata_data[i].data(),
                  field_size * sizeof(double));
    }
  }

  // Create the struct
  C::vtkPolyData c_str = { .points_state = C::SEALED,
                           .points = points,
                           .n_points = number_of_points(),
                           .m_points = number_of_points(),
                           .connectivity_state = C::SEALED,
                           .vertices_connectivity = vertices_connectivity,
                           .m_vertices_connectivity = n_vertices_connectivity,
                           .n_vertices_connectivity = n_vertices_connectivity,
                           .vertices_offsets = vertices_offsets,
                           .n_vertices_offsets = n_vertices_offsets,
                           .m_vertices_offsets = n_vertices_offsets,
                           .lines_connectivity = lines_connectivity,
                           .n_lines_connectivity = n_lines_connectivity,
                           .m_lines_connectivity = n_lines_connectivity,
                           .lines_offsets = lines_offsets,
                           .n_lines_offsets = n_lines_offsets,
                           .m_lines_offsets = n_lines_offsets,
                           .strips_connectivity = strips_connectivity,
                           .n_strips_connectivity = n_strips_connectivity,
                           .m_strips_connectivity = n_strips_connectivity,
                           .strips_offsets = strips_offsets,
                           .n_strips_offsets = n_strips_offsets,
                           .m_strips_offsets = n_strips_offsets,
                           .polygons_connectivity = polygons_connectivity,
                           .n_polygons_connectivity = n_polygons_connectivity,
                           .m_polygons_connectivity = n_polygons_connectivity,
                           .polygons_offsets = polygons_offsets,
                           .n_polygons_offsets = n_polygons_offsets,
                           .m_polygons_offsets = n_polygons_offsets,
                           .fields_state = C::SEALED,
                           .n_pointdata = n_pointdata,
                           .m_pointdata = n_pointdata,
                           .pointdata_names = pointdata_names,
                           .pointdata_ncomp = pointdata_ncomp,
                           .pointdata_data = pointdata_data,
                           .n_celldata = n_celldata,
                           .m_celldata = n_celldata,
                           .celldata_names = celldata_names,
                           .celldata_ncomp = celldata_ncomp,
                           .celldata_data = celldata_data };

  return c_str;
}

bool
vtkPolyData::points_is_sealed()
{
  return points_state == State::SEALED;
}

bool
vtkPolyData::connectivity_is_sealed()
{
  return connectivity_state == State::SEALED;
}

bool
vtkPolyData::fields_is_sealed()
{
  return fields_state == State::SEALED;
}

bool
vtkPolyData::point_exists(int64_t point)
{
  return 0 <= point && point < number_of_points();
}

void
vtkPolyData::on_add_points()
{
  if (points_is_sealed())
    ELFF_ABORT("New points may not be added after adding connectivity.\n");
}

void
vtkPolyData::on_add_connectivity()
{
  points_state = State::SEALED;
  if (connectivity_is_sealed())
    ELFF_ABORT(
      "New connectivity may not be added after adding point or cell data.\n");
}

void
vtkPolyData::on_add_field_data()
{
  points_state = State::SEALED;
  connectivity_state = State::SEALED;
  if (fields_is_sealed())
    ELFF_ABORT("New fields may not be added.\n");
}

const size_t
vtkPolyData::number_of_points() const
{
  return points.size() / 3;
}

const size_t
vtkPolyData::number_of_vertices() const
{
  return vertices_offsets.size() - 1;
}

const size_t
vtkPolyData::number_of_lines() const
{
  return lines_offsets.size() - 1;
}

const size_t
vtkPolyData::number_of_strips() const
{
  return strips_offsets.size() - 1;
}

const size_t
vtkPolyData::number_of_polygons() const
{
  return polygons_offsets.size() - 1;
}

const size_t
vtkPolyData::number_of_cells() const
{
  return number_of_vertices() + number_of_lines() + number_of_strips() +
         number_of_polygons();
}

const size_t
vtkPolyData::number_of_pointdata() const
{
  return pointdata_data.size();
}

const size_t
vtkPolyData::number_of_celldata() const
{
  return celldata_data.size();
}

void
vtkPolyData::reserve_points(size_t n)
{
  points.reserve(n * 3);
}

void
vtkPolyData::reserve_vertices(size_t n)
{
  vertices_connectivity.reserve(n);
  vertices_offsets.reserve(n + 1);
}

void
vtkPolyData::reserve_lines(size_t n)
{
  lines_connectivity.reserve(n);
  lines_offsets.reserve(n + 1);
}

void
vtkPolyData::reserve_strips(size_t n)
{
  strips_connectivity.reserve(n);
  strips_offsets.reserve(n + 1);
}

void
vtkPolyData::reserve_polygons(size_t n)
{
  polygons_connectivity.reserve(n);
  polygons_offsets.reserve(n + 1);
}

int64_t
vtkPolyData::add_pointdata_scalar(const std::string& name)
{
  return add_pointdata_vector(name, 1);
}

int64_t
vtkPolyData::add_pointdata_vector(const std::string& name, size_t ncomp)
{
  on_add_field_data();

  if (ncomp == 0) {
    ELFF_ABORT("Point-data vector must have at least one component.\n");
  }

  pointdata_names.push_back(name);
  pointdata_ncomp.push_back(ncomp);
  pointdata_data.emplace_back(number_of_points() * ncomp, 0.);

  return static_cast<int64_t>(pointdata_data.size() - 1);
}

std::vector<double>&
vtkPolyData::get_pointdata(int64_t field)
{
  const size_t id = static_cast<size_t>(field);

  if (id >= pointdata_data.size()) {
    ELFF_ABORT("Point-data field does not exist.\n");
  }

  return pointdata_data[id];
}

const std::vector<double>&
vtkPolyData::get_pointdata(int64_t field) const
{
  const size_t id = static_cast<size_t>(field);

  if (id >= pointdata_data.size()) {
    ELFF_ABORT("Point-data field does not exist.\n");
  }

  return pointdata_data[id];
}

void
vtkPolyData::set_pointdata_vector3(int64_t field,
                                   size_t point_id,
                                   const std::array<double, 3>& value)
{
  const size_t id = static_cast<size_t>(field);

  if (id >= pointdata_data.size()) {
    ELFF_ABORT("Point-data field does not exist.\n");
  }

  if (pointdata_ncomp[id] != 3) {
    ELFF_ABORT("Point-data field is not 3-component.\n");
  }

  if (point_id >= number_of_points()) {
    ELFF_ABORT("Point-data point index out of range.\n");
  }

  auto& data = pointdata_data[id];
  data[(point_id * 3) + 0] = value[0];
  data[(point_id * 3) + 1] = value[1];
  data[(point_id * 3) + 2] = value[2];
}

int64_t
vtkPolyData::add_celldata_scalar(const std::string& name)
{
  return add_celldata_vector(name, 1);
}

int64_t
vtkPolyData::add_celldata_vector(const std::string& name, size_t ncomp)
{
  on_add_field_data();

  if (ncomp == 0) {
    ELFF_ABORT("Cell-data vector must have at least one component.\n");
  }

  celldata_names.push_back(name);
  celldata_ncomp.push_back(ncomp);
  celldata_data.emplace_back(number_of_cells() * ncomp, 0.);

  return static_cast<int64_t>(celldata_data.size() - 1);
}

std::vector<double>&
vtkPolyData::get_celldata(int64_t field)
{
  const size_t id = static_cast<size_t>(field);

  if (id >= celldata_data.size()) {
    ELFF_ABORT("Cell-data field does not exist.\n");
  }

  return celldata_data[id];
}

const std::vector<double>&
vtkPolyData::get_celldata(int64_t field) const
{
  const size_t id = static_cast<size_t>(field);

  if (id >= celldata_data.size()) {
    ELFF_ABORT("Cell-data field does not exist.\n");
  }

  return celldata_data[id];
}

void
vtkPolyData::set_celldata_vector3(int64_t field,
                                  size_t cell_id,
                                  const std::array<double, 3>& value)
{
  const size_t id = static_cast<size_t>(field);

  if (id >= celldata_data.size()) {
    ELFF_ABORT("Cell-data field does not exist.\n");
  }

  if (celldata_ncomp[id] != 3) {
    ELFF_ABORT("Cell-data field is not 3-component.\n");
  }

  if (cell_id >= number_of_cells()) {
    ELFF_ABORT("Cell-data cell index out of range.\n");
  }

  auto& data = celldata_data[id];
  data[(cell_id * 3) + 0] = value[0];
  data[(cell_id * 3) + 1] = value[1];
  data[(cell_id * 3) + 2] = value[2];
}

int64_t
vtkPolyData::add_point(float x, float y, float z)
{
  on_add_points();

  points.push_back(x);
  points.push_back(y);
  points.push_back(z);

  return number_of_points() - 1;
}

int64_t
vtkPolyData::add_vertex(int64_t vertex_point)
{
  on_add_connectivity();

  if (!point_exists(vertex_point)) {
    ELFF_ABORT("Point does not exist.\n");
  }

  vertices_connectivity.push_back(vertex_point);
  vertices_offsets.push_back(vertices_connectivity.size());

  return number_of_vertices() - 1;
}

int64_t
vtkPolyData::add_poly_vertex(std::vector<int64_t> poly_vertex_points)
{
  on_add_connectivity();

  if (poly_vertex_points.size() < 2) {
    ELFF_ABORT("Poly vertex must contain at least 2 points.\n");
  }

  for (auto point : poly_vertex_points) {
    if (!point_exists(point)) {
      ELFF_ABORT("One of the points does not exist.\n");
    }
  }

  vertices_connectivity.insert(vertices_connectivity.end(),
                               poly_vertex_points.begin(),
                               poly_vertex_points.end());
  vertices_offsets.push_back(vertices_connectivity.size());

  return number_of_vertices() - 1;
}

int64_t
vtkPolyData::add_line(int64_t point_1, int64_t point_2)
{
  on_add_connectivity();

  if (!point_exists(point_1) || !point_exists(point_2)) {
    ELFF_ABORT("One of the points does not exist.\n");
  }

  lines_connectivity.push_back(point_1);
  lines_connectivity.push_back(point_2);
  lines_offsets.push_back(lines_connectivity.size());

  return number_of_lines() - 1;
}

int64_t
vtkPolyData::add_polyline(std::vector<int64_t> polyline_points)
{
  on_add_connectivity();

  if (polyline_points.size() < 2) {
    ELFF_ABORT("Line must contain at least 2 points.\n");
  }

  for (auto point : polyline_points) {
    if (!point_exists(point)) {
      ELFF_ABORT("One of the points does not exist.\n");
    }
  }

  lines_connectivity.insert(
    lines_connectivity.end(), polyline_points.begin(), polyline_points.end());
  lines_offsets.push_back(lines_connectivity.size());

  return number_of_lines() - 1;
}

int64_t
vtkPolyData::add_polygon(std::vector<int64_t> polygon_points)
{
  on_add_connectivity();

  if (polygon_points.size() < 3) {
    ELFF_ABORT("Polygon must contain at least 3 points.\n");
  }

  for (auto point : polygon_points) {
    if (!point_exists(point)) {
      ELFF_ABORT("One of the points does not exist.\n");
    }
  }

  polygons_connectivity.insert(
    polygons_connectivity.end(), polygon_points.begin(), polygon_points.end());
  polygons_offsets.push_back(polygons_connectivity.size());

  return number_of_polygons() - 1;
}

void
vtkPolyData::append(const vtkPolyData& other)
{
  if (this == &other) {
    vtkPolyData copy = other;
    append(copy);
    return;
  }

  ELFF_ASSERT(pointdata_names.size() == other.pointdata_names.size(),
              "vtkPolyData::append(): point-data field count mismatch.\n");
  for (size_t i = 0; i < pointdata_names.size(); ++i) {
    ELFF_ASSERT(pointdata_names[i] == other.pointdata_names[i],
                "vtkPolyData::append(): point-data field name mismatch.\n");
    ELFF_ASSERT(pointdata_ncomp[i] == other.pointdata_ncomp[i],
                "vtkPolyData::append(): point-data component count mismatch.\n");
  }

  ELFF_ASSERT(celldata_names.size() == other.celldata_names.size(),
              "vtkPolyData::append(): cell-data field count mismatch.\n");
  for (size_t i = 0; i < celldata_names.size(); ++i) {
    ELFF_ASSERT(celldata_names[i] == other.celldata_names[i],
                "vtkPolyData::append(): cell-data field name mismatch.\n");
    ELFF_ASSERT(celldata_ncomp[i] == other.celldata_ncomp[i],
                "vtkPolyData::append(): cell-data component count mismatch.\n");
  }

  C::vtkPolyData c_dst = this->to_c_struct();
  C::vtkPolyData c_src = other.to_c_struct();

  C::vtk_polydata_append(&c_dst, &c_src);

  points_state = (c_dst.points_state == C::SEALED) ? State::SEALED : State::BUILDING;
  points.assign(c_dst.points, c_dst.points + (c_dst.n_points * 3));

  connectivity_state =
    (c_dst.connectivity_state == C::SEALED) ? State::SEALED : State::BUILDING;

  vertices_connectivity.assign(c_dst.vertices_connectivity,
                               c_dst.vertices_connectivity + c_dst.n_vertices_connectivity);
  vertices_offsets.assign(
    c_dst.vertices_offsets, c_dst.vertices_offsets + c_dst.n_vertices_offsets);

  lines_connectivity.assign(
    c_dst.lines_connectivity, c_dst.lines_connectivity + c_dst.n_lines_connectivity);
  lines_offsets.assign(c_dst.lines_offsets, c_dst.lines_offsets + c_dst.n_lines_offsets);

  strips_connectivity.assign(
    c_dst.strips_connectivity, c_dst.strips_connectivity + c_dst.n_strips_connectivity);
  strips_offsets.assign(c_dst.strips_offsets, c_dst.strips_offsets + c_dst.n_strips_offsets);

  polygons_connectivity.assign(c_dst.polygons_connectivity,
                               c_dst.polygons_connectivity + c_dst.n_polygons_connectivity);
  polygons_offsets.assign(
    c_dst.polygons_offsets, c_dst.polygons_offsets + c_dst.n_polygons_offsets);

  fields_state = (c_dst.fields_state == C::SEALED) ? State::SEALED : State::BUILDING;

  pointdata_names.resize(c_dst.n_pointdata);
  pointdata_ncomp.resize(c_dst.n_pointdata);
  pointdata_data.resize(c_dst.n_pointdata);

  for (size_t i = 0; i < c_dst.n_pointdata; ++i) {
    pointdata_names[i] = c_dst.pointdata_names[i];
    pointdata_ncomp[i] = c_dst.pointdata_ncomp[i];
    pointdata_data[i].assign(c_dst.pointdata_data[i],
                             c_dst.pointdata_data[i] +
                               (c_dst.n_points * c_dst.pointdata_ncomp[i]));
  }

  const size_t n_cells =
    (c_dst.n_vertices_offsets - 1) + (c_dst.n_lines_offsets - 1) +
    (c_dst.n_strips_offsets - 1) + (c_dst.n_polygons_offsets - 1);

  celldata_names.resize(c_dst.n_celldata);
  celldata_ncomp.resize(c_dst.n_celldata);
  celldata_data.resize(c_dst.n_celldata);

  for (size_t i = 0; i < c_dst.n_celldata; ++i) {
    celldata_names[i] = c_dst.celldata_names[i];
    celldata_ncomp[i] = c_dst.celldata_ncomp[i];
    celldata_data[i].assign(c_dst.celldata_data[i],
                            c_dst.celldata_data[i] +
                              (n_cells * c_dst.celldata_ncomp[i]));
  }

  C::vtk_polydata_free(&c_dst);
  C::vtk_polydata_free(&c_src);
}

vtkPolyData
vtkPolyData::append_many(const std::vector<vtkPolyData>& datasets)
{
  vtkPolyData out;

  if (datasets.empty()) {
    return out;
  }

  out = datasets.front();
  for (size_t i = 1; i < datasets.size(); ++i) {
    out.append(datasets[i]);
  }

  return out;
}

}
}
}
