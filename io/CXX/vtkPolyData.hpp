#pragma once

#include "config/config.hpp"
#include "general/error.hpp"

#include "io/C/vtkPolyData.hpp"

#include <vector>

namespace ELFF {
namespace io {
namespace CXX {

class vtkHDFPolyData;

/**
 * @brief Wrapper for the @ref C::vtkPolyData class
 *
 * This class actually implements its own logic to ensure correct vtkPolyData
 * fields, but is easily convertible to the C class.
 */
class vtkPolyData
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

  /**
   * @brief Returns true if the points are sealed
   */
  bool points_is_sealed();

  /**
   * @brief Returns true if the connectivities are sealed
   */
  bool connectivity_is_sealed();

  /**
   * @brief Returns true if the fields are sealed
   */
  bool fields_is_sealed();

  /**
   * @brief Returns true if the point exists
   */
  bool point_exists(int64_t point);

  /**
   * @brief Called on addition of points
   */
  void on_add_points();

  /**
   * @brief Called on addition of connectivity data
   */
  void on_add_connectivity();

  /**
   * @brief Called on addition of field data
   */
  void on_add_field_data();

  /**
   * @brief Return a C vtkPolyData struct with heap arrays belonging to it
   */
  C::vtkPolyData to_c_struct();

  friend class vtkHDFPolyData;

public:
  vtkPolyData()
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
    , fields_state(State::BUILDING) {};

  /**
   * @brief Returns the number of points
   */
  const size_t number_of_points() const;

  /**
   * @brief Returns the number of vertices
   */
  const size_t number_of_vertices() const;

  /**
   * @brief Returns the number of lines
   */
  const size_t number_of_lines() const;

  /**
   * @brief Returns the number of strips
   */
  const size_t number_of_strips() const;

  /**
   * @brief Returns the number of polygons
   */
  const size_t number_of_polygons() const;

  /**
   * @brief Pre-allocation for points data
   */
  void reserve_points(size_t n);

  /**
   * @brief Pre-allocation for vertices cells
   */
  void reserve_vertices(size_t n);

  /**
   * @brief Pre-allocation for lines cells
   */
  void reserve_lines(size_t n);

  /**
   * @brief Pre-allocation for strips cells
   */
  void reserve_strips(size_t n);

  /**
   * @brief Pre-allocation for polygons
   */
  void reserve_polygons(size_t n);

  /**
   * @brief Add a new point to the dataset
   *
   * @param x
   * @param y
   * @param z
   */
  int64_t add_point(float x, float y, float z);

  /**
   * @brief Add a new vertex cell
   *
   * @param vertex_point The point index for the new vertex
   */
  int64_t add_vertex(int64_t vertex_point);

  /**
   * @brief Add a poly vertex cell
   *
   * @param poly_vertex_points The point index for the new vertex
   */
  int64_t add_poly_vertex(std::vector<int64_t> poly_vertex_points);

  /**
   * @brief Add a line cell
   *
   * @param point_1
   * @param point_2
   */
  int64_t add_line(int64_t point_1, int64_t point_2);

  /**
   * @brief Add a polyline cell
   *
   * @param polyline_points
   */
  int64_t add_polyline(std::vector<int64_t> polyline_points);

  /**
   * @brief Add a polygon cell
   *
   * @param polygon_points
   */
  int64_t add_polygon(std::vector<int64_t> polygon_points);
};

}
}
}