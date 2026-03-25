#pragma once

#include "elff/config/config.hpp"
#include "elff/models/rope/RopeMesh.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

class Rope
{
public:
  typedef enum
  {
    free_bc,
    simple_bc,
    point_force_bc,
  } RopeBCType;

  typedef enum
  {
    left,
    right,
  } RopeBCEnd;

  typedef struct
  {
    double position[3];
    double slope[3];
    double force[3];
    double torque[3];
  } RopeBCVals;

  typedef struct
  {
    RopeBCEnd end[2];
    RopeBCType type[2];
    RopeBCVals vals[2];
  } RopeBCs;

  typedef enum
  {
    INVALID = -1,
    NONE,
    UNIFORM,
    NONUNIFORM,
  } RopeLoadType;

  virtual void solve();
  virtual void solve(std::array<real_t, 3> uniform_load);
  virtual void solve(std::vector<std::array<real_t, 3>> nonuniform_load);
  virtual void solve(real_t dt);
  virtual void solve(real_t dt, std::array<real_t, 3> uniform_load);
  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> nonuniform_load);
  virtual void apply_initial_condition(RopeMesh& mesh) = 0;
  virtual RopeMesh& get_mesh();
  virtual const RopeMesh& get_mesh() const;

protected:
  bool is_time_dependent;
  size_t time_iter;
  real_t t;
  size_t dimension;
  real_t mu;
  RopeMesh mesh;
  RopeBCs boundary_conditions;
  Rope();
  Rope(real_t length, size_t nodes, RopeBCs bcs);
  Rope(real_t length, real_t mu, size_t nodes, RopeBCs bcs);
};

} // namespace Models
} // namespace ELFF
