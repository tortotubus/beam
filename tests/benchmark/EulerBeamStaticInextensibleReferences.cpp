#include "EulerBeamStaticInextensibleReferences.hpp"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

namespace ELFF {

TEST(EulerBeamStaticInextensibleReferences, BisshoppAndDrucker)
{
  const std::string csv_filename = "bisshopp-drucker.csv";
  const std::string gp_filename = "bisshopp-drucker.gp";
  const std::string tex_filename = "bisshopp-drucker.tex";

  std::ofstream csv(csv_filename);
  std::ofstream gp(gp_filename);

  ASSERT_TRUE(csv.is_open());
  ASSERT_TRUE(gp.is_open());

  const double EI = 1.;
  const double length = 1.;
  const double p_min = 0.;
  const double p_max = 10.;
  const double p_delta = 0.1;

  double p_i = p_min;

  std::vector<double> nondim_x_tip;
  std::vector<double> nondim_y_displacement;
  std::vector<double> nondim_load;

  while (p_i <= p_max) {
    const auto res_i = BisshoppAndDrucker1945(length, EI, p_i);
    const double nondim_x_tip_i = (length - res_i.A) / length;
    const double nondim_y_disp_i = res_i.delta / length;
    const double nondim_load_i = (p_i * length * length) / EI;

    nondim_x_tip.push_back(nondim_x_tip_i);
    nondim_y_displacement.push_back(nondim_y_disp_i);
    nondim_load.push_back(nondim_load_i);

    p_i += p_delta;
  }

  csv << "x_tip,y_disp,load\n";
  for (size_t i = 0; i < nondim_load.size(); ++i) {
    csv << nondim_x_tip[i] << "," << nondim_y_displacement[i] << ","
        << nondim_load[i] << "\n";
  }

  gp << "if (!exists(\"csv_file\")) csv_file = \"" << csv_filename << "\"\n"
     << "if (!exists(\"plot_output\")) plot_output = \"../vector/"
     << tex_filename << "\"\n";

  gp << R"gnuplot(
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel ""
set ylabel "$PL^2/B$"
set xtics 0.0,0.1,1.00
set ytics 0.0,1.0,10.0
set xrange [0:1]
set yrange [0:10]
set grid
set size square
set key outside right 
plot \
  csv_file every ::1 using 1:3 with lines title "$\\frac{x_{\\rm tip}}{L}$" smooth csplines, \
  csv_file every ::1 using 2:3 with lines title "$\\frac{y_{\\rm tip}}{L}$" smooth csplines
)gnuplot";
}

} // namespace ELFF
