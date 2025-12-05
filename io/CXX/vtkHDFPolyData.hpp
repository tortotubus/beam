#pragma once

#include "io/C/vtkPolyData.hpp"

#include <string>

namespace beam {
namespace io {
namespace CXX {

class vtkPolyData;

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

  void write_new_static(bool overwrite);
  void write_new_static(bool overwrite, int compression_level);

  void write_new_transient(bool overwrite, float time);
  void write_new_transient(bool overwrite, float time, int compression_level);

  void append_transient(float time);
  void append_transient(float time, int compression_level);
};

}
}
}