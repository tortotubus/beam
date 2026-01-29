#include "EulerBeamStaticInextensibleReferences.hpp"

#include <gtest/gtest.h>

namespace ELFF {

TEST(EulerBeamStaticInextensibleReferences, BisshoppAndDrucker) {
  double EI = 1., length = 1.;
  double p_min = 0., p_max = 10., p_delta = 0.1;

  double p_i = p_min;

  std::vector<double> nondim_x_displacement;
  std::vector<double> nondim_y_displacement;
  std::vector<double> nondim_load;

  while (p_i <= p_max) {
    auto res_i = BisshoppAndDrucker1945(length, EI, p_i);
    double nondim_y_disp_i = res_i.delta / length;
    double nondim_x_disp_i = (length - res_i.A) / length;
    nondim_x_displacement.push_back(nondim_x_disp_i);
    nondim_y_displacement.push_back(nondim_y_disp_i);
    nondim_load.push_back((p_i*length*length)/EI);
    p_i += p_delta;
  }

  
  FILE *pipe = popen("gnuplot -persist", "w");

  if (!pipe) {
    return;
  }

  // Configure the plot
  fprintf(pipe, "set title 'Bisshopp and Drucker 1945' \n");
  fprintf(pipe, "set terminal qt enhanced \n");
  fprintf(pipe, "set xlabel '' \n");
  fprintf(pipe, "set ylabel 'PL^2/(EI)' \n");
  fprintf(pipe, "set xtics 0.0, 0.1, 1.00\n");
  fprintf(pipe, "set ytics 0.0, 1.0, 10.0\n");
  // fprintf(pipe, "set xrange [%f:%f]\n", 10., 0.);
  // fprintf(pipe, "set yrange [%f:%f]\n", 10., 0.);
  fprintf(pipe, "set grid\n");
  fprintf(pipe, "set size square\n");      // <= disable extension to tics :contentReference[oaicite:0]{index=0}  
  fprintf(pipe, "plot '-' using 1:2 with lines title '(L-X_{tip})/L' smooth csplines, ");
  fprintf(pipe, "'-' using 1:2 with lines title 'y_{tip}/L' smooth csplines\n");

  // Send the data points
  for (size_t i = 0; i < nondim_x_displacement.size(); ++i) {
    fprintf(pipe, "%f %f\n", nondim_x_displacement[i], nondim_load[i]);
  }
  fprintf(pipe, "e\n"); // End of data marker for gnuplot

  // Send the data points
  for (size_t i = 0; i < nondim_y_displacement.size(); ++i) {
    fprintf(pipe, "%f %f\n", nondim_y_displacement[i], nondim_load[i]);
  }
  fprintf(pipe, "e\n"); // End of data marker for gnuplot

  // Clean up
  pclose(pipe);
  
};

}