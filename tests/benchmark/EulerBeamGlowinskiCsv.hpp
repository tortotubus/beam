#pragma once

#include <elff/models/beam/EulerBeamMesh.hpp>

#include <iomanip>
#include <limits>
#include <ostream>

namespace ELFF {
namespace Benchmark {

inline void
write_glowinski_csv_header(std::ostream& os)
{
  os << "# time,step,node,s,x,y,z,vx,vy,vz\n";
  os << std::setprecision(std::numeric_limits<real_t>::max_digits10);
}

inline void
write_glowinski_csv_frame(std::ostream& os,
                          size_t step,
                          real_t time,
                          Models::EulerBeamMesh& mesh)
{
  auto& centerline = mesh.get_centerline();
  auto& velocity = mesh.get_centerline_velocity();

  for (size_t node = 0; node < mesh.get_nodes(); ++node) {
    os << time << ',' << step << ',' << node << ','
       << mesh.get_curvilinear_axis(node) << ',' << centerline[node][0] << ','
       << centerline[node][1] << ',' << centerline[node][2] << ','
       << velocity[node][0] << ',' << velocity[node][1] << ','
       << velocity[node][2] << '\n';
  }

  os << '\n';
}

} // namespace Benchmark
} // namespace ELFF
