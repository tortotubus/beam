#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <gtest/gtest.h>

#include "EulerBeamDynamicInextensibleReferences.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamDynamicInextensibleReferences, MMF1)
{
  const real_t dt = 0.02;
  const real_t tf = 1.;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const int N = 240;

  for (int ti = -2; ti < Nt; ti++) {

    const std::string filename = "mmf1.vtkhdf";
    double t = ti * static_cast<double>(dt);
    auto mesh = ManufacturedDynamicResult1(N, t);
    vtkPolyData pd = mesh.to_vtk_polydata();
    vtkHDFPolyData hdf_pd(filename, pd);

    if (ti == 0) {
      hdf_pd.write_new_transient(true, t);
    } else {
      hdf_pd.append_transient(t);
    }
  }
}

};