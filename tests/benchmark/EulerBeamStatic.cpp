#include <gtest/gtest.h>
#include "models/beam/EulerBeamStatic.hpp"

using beam::EulerBeam;
using beam::EulerBeamBCs;
using beam::EulerBeamMesh;
using beam::EulerBeamStatic;
using beam::real_t;
using beam::left;
using beam::right;
using beam::clamped_bc;
using beam::free_bc;
using beam::point_force_bc;
using beam::point_torque_bc;


