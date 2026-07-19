#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <gtest/gtest.h>

#include <fstream>
#include <string>

#include "EulerBeamDynamicInextensibleReferences.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamDynamicInextensibleReferences, MMF1)
{
  const std::string hdf_filename = "mmf1.vtkhdf";
  const std::string csv_filename = "mmf1.csv";
  const std::string gp_filename = "mmf1.gp";
  const std::string tex_filename = "mmf1.tex";

  std::ofstream csv(csv_filename);
  std::ofstream gp(gp_filename);

  ASSERT_TRUE(csv.is_open());
  ASSERT_TRUE(gp.is_open());

  const real_t dt = 0.02;
  const real_t tf = 1.;
  const int Nt = static_cast<int>(std::ceil(tf / dt));
  const int N = 240;

  csv << "# t,s,x,y,z,sx,sy,sz,vx,vy,vz,ax,ay,az\n";

  gp << "if (!exists(\"csv_file\")) csv_file = \"" << csv_filename << "\"\n"
     << "if (!exists(\"plot_output\")) plot_output = \"../vector/"
     << tex_filename << "\"\n";

  gp << R"gnuplot(
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$x$"
set ylabel "$y$"
set grid
set size square
set key outside right
plot \
  csv_file index 0 using 3:4 with lines title "$t=-0.04$", \
  csv_file index 2 using 3:4 with lines title "$t=0$", \
  csv_file index 27 using 3:4 with lines title "$t=0.5$", \
  csv_file index 51 using 3:4 with lines title "$t=0.98$"
)gnuplot";

  for (int ti = -2; ti < Nt; ti++) {
    double t = ti * static_cast<double>(dt);
    auto mesh = ManufacturedDynamicResult1(N, t);

    const auto& s = mesh.get_curvilinear_axis();
    const auto& position = mesh.get_centerline();
    const auto& slope = mesh.get_slope();
    const auto& velocity = mesh.get_centerline_velocity();
    const auto& acceleration = mesh.get_centerline_acceleration();

    for (std::size_t i = 0; i < mesh.get_nodes(); ++i) {
      csv << t << "," << s[i] << "," << position[i][0] << ","
          << position[i][1] << "," << position[i][2] << "," << slope[i][0]
          << "," << slope[i][1] << "," << slope[i][2] << ","
          << velocity[i][0] << "," << velocity[i][1] << ","
          << velocity[i][2] << "," << acceleration[i][0] << ","
          << acceleration[i][1] << "," << acceleration[i][2] << "\n";
    }
    csv << "\n\n";

    vtkPolyData pd = mesh.to_vtk_polydata();
    vtkHDFPolyData hdf_pd(hdf_filename, pd);

    if (ti == -2) {
      hdf_pd.write_new_transient(true, t);
    } else {
      hdf_pd.append_transient(t);
    }
  }
}

};
