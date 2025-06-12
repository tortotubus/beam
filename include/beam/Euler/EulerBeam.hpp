#pragma once

namespace beam {

typedef enum {
  free_bc,
  simple_bc,
  clamped_bc,
} EulerBeamBCType;

typedef enum {
  left,
  right,
} EulerBeamBCEnd;

typedef struct {
  double position[3];
  double slope[3];
} EulerBeamBCVals;

typedef struct {
  EulerBeamBCEnd end[2];
  EulerBeamBCType type[2];
  EulerBeamBCVals vals[2];
} EulerBeamBCs;

} // namespace beam
