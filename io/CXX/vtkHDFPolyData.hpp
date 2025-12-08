#pragma once

#include "io/C/vtkPolyData.hpp"

#include <string>

namespace beam {
namespace io {
namespace CXX {


class vtkPolyData;

/**
 * @brief Wrapper for the @ref C::vtkHDFPolyData class
 */
class vtkHDFPolyData
{
private:
  std::string filename;
  C::vtkPolyData vtk_polydata;

public:
  vtkHDFPolyData(std::string filename);
  vtkHDFPolyData(std::string filename, C::vtkPolyData c_vtk_polydata);
  vtkHDFPolyData(std::string filename, CXX::vtkPolyData& cxx_vtk_polydata);

  ~vtkHDFPolyData();

  /**
   * @brief Write a new time-independent file
   * 
   * @param overwrite If set to true, any existing files will be overwritten
   */
  void write_new_static(bool overwrite);

  /**
   * @brief Write a new time-independent file
   * 
   * @param overwrite If set to true, any existing files will be overwritten
   * @param compression_level Level of compression [0,9] to be used
   */
  void write_new_static(bool overwrite, int compression_level);

  /**
   * @brief Write a new time-dependent file
   * 
   * @param overwrite If set to true, any existing files will be overwritten
   * @param time The time value to write 
   */
  void write_new_transient(bool overwrite, float time);

  /**
   * @brief Write a new time-dependent file
   * 
   * @param overwrite If set to true, any existing files will be overwritten
   * @param time The time value to write 
   * @param compression_level Level of compression [0,9] to be used
   */
  void write_new_transient(bool overwrite, float time, int compression_level);

  /**
   * @brief Append to an existing time-dependent file
   * 
   * @param time The time value to write
   */
  void append_transient(float time);

  /**
   * @brief Append to an existing time-dependent file
   * 
   * @param time The time value to write
   * @param compression_level Level of compression [0,9] to be used
   */
  void append_transient(float time, int compression_level);
};

}
}
}