#include "io/CXX/vtkPolyData.hpp"

#include <cstdlib>   // malloc, free
#include <cstring>   // std::memcpy

namespace ELFF {
namespace io {
namespace CXX {

C::vtkPolyData
vtkPolyData::to_c_struct()
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
                           .fields_state = C::SEALED };

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

}
}
}