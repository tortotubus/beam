
#include <gtest/gtest.h>
#include "models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

#include "EulerBeamStaticInextensibleReferences.hpp"

namespace ELFF {
using namespace io::CXX;

TEST(EulerBeamStaticInextensibleMoMSparseTest, BisshoppAndDrucker) {

  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e5;
  size_t nodes = 40;

  real_t tip_force_y = -1;

  double comparison_tol = 1e-5;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = {EulerBeam::left, EulerBeam::right},
    .type = {EulerBeam::clamped_bc, EulerBeam::point_force_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {
      .force = {0,tip_force_y,0}
    }}
  };

  EulerBeamStaticInextensibleMoMSparse static_beam(length, EI, nodes, boundary_conditions, r_pentalty);

  static_beam.solve();
  EulerBeamMesh& mesh = static_beam.get_mesh(); 
  auto centerline = mesh.get_centerline();
  std::array<real_t, 3> tip = centerline[nodes-1];

  BisshoppAndDrucker1945Result res = BisshoppAndDrucker1945(length, EI, -tip_force_y);

  EXPECT_NEAR(std::abs(length-tip[0]), res.A, comparison_tol);
  EXPECT_NEAR(std::abs(tip[1]), res.delta, comparison_tol);

  FILE *pipe = popen("gnuplot -persist", "w");
  if (!pipe) {
    ELFF_ABORT("Failed to open pipe to gnuplot");
  }

  // Configure the plot
  fprintf(pipe, "set title 'Bisshopp and Drucker Comparison'\n");
  fprintf(pipe, "set xlabel 'x'\n");
  fprintf(pipe, "set ylabel 'y'\n");
  fprintf(pipe, "set xrange [%f:%f]\n",         -.1, length + .1);
  fprintf(pipe, "set yrange [%f:%f]\n", -length -.1,          .1);
  fprintf(pipe, "set grid\n");
  fprintf(pipe, "set size square\n");      // <= disable extension to tics :contentReference[oaicite:0]{index=0}  
  fprintf(pipe, "plot '-' using 1:2 with lines title 'beam' smooth csplines, ");
  fprintf(pipe, "'-' using 1:2 with points pt 2 ps 2 lc rgb 'red' title 'ref'\n");
  // Send the data points
  for (size_t i = 0; i < centerline.size(); ++i) {
    fprintf(pipe, "%f %f\n", centerline[i][0], centerline[i][1]);
  }
  fprintf(pipe, "e\n"); // End of data marker for gnuplot

  // Send the reference point (second dataset) -> red X at (res.A, res.delta)
  fprintf(pipe, "%f %f\n", length-res.A, -res.delta);
  fprintf(pipe, "e\n"); // end second dataset

  // Clean up
  pclose(pipe);

};

}