#pragma once

#include "config/config.hpp"
#include "general/error.hpp"

#include <map>
#include <string>
#include <vector>

namespace beam {

/**
 *
 */
class VTKPolyData
{
protected:
  enum class State
  {
    SEALED,
    BUILDING,
  };

  State points_state;
  std::vector<float> points;

  State connectivity_state;

  std::vector<int64_t> vertices_connectivity;
  std::vector<int64_t> vertices_offsets;

  std::vector<int64_t> lines_connectivity;
  std::vector<int64_t> lines_offsets;

  std::vector<int64_t> strips_connectivity;
  std::vector<int64_t> strips_offsets;

  std::vector<int64_t> polygons_connectivity;
  std::vector<int64_t> polygons_offsets;

  State fields_state;
  std::map<std::string, size_t> fields_name_map;
  // std::vector<PolyDataField> fields;

  bool points_is_sealed() { return points_state == State::SEALED; }
  bool connectivity_is_sealed() { return connectivity_state == State::SEALED; }
  bool fields_is_sealed() { return fields_state == State::SEALED; }

  bool point_exists(int64_t point) { return 0 <= point && point < number_of_points(); }

  void on_add_points()
  {
    if (points_is_sealed())
      BEAM_ABORT("New points may not be added after adding connectivity.\n");
  }

  void on_add_connectivity()
  {
    points_state = State::SEALED;
    if (connectivity_is_sealed())
      BEAM_ABORT(
        "New connectivity may not be added after adding point or cell data.\n");
  }

  void on_add_field_data()
  {
    points_state = State::SEALED;
    connectivity_state = State::SEALED;
    if (fields_is_sealed())
      BEAM_ABORT("New fields may not be added.\n");
  }

public:
  VTKPolyData()
    : points_state(State::BUILDING)
    , points(0)
    , connectivity_state(State::BUILDING)
    , vertices_connectivity(0)
    , vertices_offsets({ 0 })
    , lines_connectivity(0)
    , lines_offsets({ 0 })
    , strips_connectivity(0)
    , strips_offsets({ 0 })
    , polygons_connectivity(0)
    , polygons_offsets({ 0 })
    , fields_state(State::BUILDING)
    , fields_name_map()
    // , fields()
  {
  }

  friend class VTKHDFPolyData;

  /**
   * 
   */
  const size_t number_of_points() const {
    return points.size() / 3;
  }

  /**
   * 
   */
  const size_t number_of_vertices() const {
    return vertices_offsets.size() - 1;
  }

  /**
   * 
   */
  const size_t number_of_lines() const {
    return lines_offsets.size() - 1;
  }

  /**
   * 
   */
  const size_t number_of_strips() const {
    return strips_offsets.size() - 1;
  }

  /**
   * 
   */
  const size_t number_of_polygons() const {
    return polygons_offsets.size() - 1;
  }

  /**
   *
   */
  void reserve_points(size_t n) { points.reserve(n * 3); }

  /**
   *
   */
  void reserve_vertices(size_t n)
  {
    vertices_connectivity.reserve(n);
    vertices_offsets.reserve(n+1);
  }

  /**
   *
   */
  void reserve_lines(size_t n)
  {
    lines_connectivity.reserve(n);
    lines_offsets.reserve(n+1);
  }

  /**
   *
   */
  void reserve_strips(size_t n)
  {
    strips_connectivity.reserve(n);
    strips_offsets.reserve(n+1);
  }

  /**
   *
   */
  void reserve_polygons(size_t n)
  {
    polygons_connectivity.reserve(n);
    polygons_offsets.reserve(n+1);
  }

  /**
   *
   */
  int64_t add_point(float x, float y, float z)
  {
    on_add_points();

    points.push_back(x);
    points.push_back(y);
    points.push_back(z);

    return number_of_points() - 1;
  }

  /**
   *
   */
  int64_t add_vertex(int64_t vertex_point)
  {
    on_add_connectivity();

    if (!point_exists(vertex_point)) {
      BEAM_ABORT("Point does not exist.\n");
    }

    vertices_connectivity.push_back(vertex_point);
    vertices_offsets.push_back(vertices_connectivity.size());

    return number_of_vertices() - 1;
  }

  int64_t add_poly_vertex(std::vector<int64_t> poly_vertex_points)
  {
    on_add_connectivity();

    if (poly_vertex_points.size() < 2) {
      BEAM_ABORT("Poly vertex must contain at least 2 points.\n");
    }

    for (auto point : poly_vertex_points) {
      if (!point_exists(point)) {
        BEAM_ABORT("One of the points does not exist.\n");
      }
    }

    vertices_connectivity.insert(vertices_connectivity.end(),
                                 poly_vertex_points.begin(),
                                 poly_vertex_points.end());
    vertices_offsets.push_back(vertices_connectivity.size());

    return number_of_vertices() - 1;
  }

  /**
   *
   */
  int64_t add_line(int64_t point_1, int64_t point_2)
  {
    on_add_connectivity();

    if (!point_exists(point_1) || !point_exists(point_2)) {
      BEAM_ABORT("One of the points does not exist.\n");
    }

    lines_connectivity.push_back(point_1);
    lines_connectivity.push_back(point_2);
    lines_offsets.push_back(lines_connectivity.size());

    return number_of_lines() - 1;
  }

  /**
   *
   */
  int64_t add_polyline(std::vector<int64_t> polyline_points)
  {
    on_add_connectivity();

    if (polyline_points.size() < 2) {
      BEAM_ABORT("Line must contain at least 2 points.\n");
    }

    for (auto point : polyline_points) {
      if (!point_exists(point)) {
        BEAM_ABORT("One of the points does not exist.\n");
      }
    }

    lines_connectivity.insert(
      lines_connectivity.end(), polyline_points.begin(), polyline_points.end());
    lines_offsets.push_back(lines_connectivity.size());

    return number_of_lines() - 1;
  }

  /**
   *
   */
  int64_t add_polygon(std::vector<int64_t> polygon_points)
  {
    on_add_connectivity();

    if (polygon_points.size() < 3) {
      BEAM_ABORT("Polygon must contain at least 3 points.\n");
    }

    for (auto point : polygon_points) {
      if (!point_exists(point)) {
        BEAM_ABORT("One of the points does not exist.\n");
      }
    }

    polygons_connectivity.insert(polygons_connectivity.end(),
                                 polygon_points.begin(),
                                 polygon_points.end());
    polygons_offsets.push_back(polygons_connectivity.size());

    return number_of_polygons() - 1;
  }
};

}