// vtkPolyDataTest.cpp
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>

// Adjust include path / name as needed:
#include "io/CXX/vtkPolyData.hpp"
#include "io/CXX/vtkHDFPolyData.hpp"

namespace ELFF {
namespace io {
namespace CXX {

TEST(vtkHDFPolyDataStaticCXX, WriteLine)
{
  size_t n_points = 5;
  double ds = 0.1;

  vtkPolyData pd;

  for (size_t i = 0; i < n_points; i++) {
    pd.add_point(0 + ds * i, 0, 0);
  }

  for (size_t i = 0; i < n_points - 1; i++) {
    pd.add_line(i,i+1);
  }

  for (size_t i = 0; i < n_points; i++) {
    pd.add_vertex(i);
  }

  std::string fname = "polydata_static_CXX.vtkhdf";
  vtkHDFPolyData hdf_pd(fname, pd);
  
  bool overwrite = true;
  hdf_pd.write_new_static(overwrite);
  
}

}
}
}