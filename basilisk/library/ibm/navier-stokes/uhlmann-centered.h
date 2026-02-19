#include "run.h"
#include "timestep.h"
#include "bcg.h"

#if EMBED
#include "library/ibm/navier-stokes/viscosity-embed-ibm.h"
#else
#include "library/ibm/navier-stokes/viscosity-ibm.h"
#endif

#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBEvents.h"
#include "library/ibm/IBAdapt.h"

scalar p[];
vector u[];
scalar pf[];
face vector uf[];

vector f_ibm[];
vector a_ibm[];

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;

double alpha_split = 0.5;
double beta_split = 0.5;

event defaults (i = 0) {
  mgp = (mgstats) {0};
  mgpf = (mgstats) {0};
  mgu = (mgstats) {0};

  p.nodump = pf.nodump = true;

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  } else if (!is_constant (alpha.x)) {
    face vector alphav = alpha;
    foreach_face () alphav.x[] = fm.x[];
  }

#if TREE
  uf.x.refine = refine_face_solenoidal;
#if EMBED 
  uf.x.refine = refine_face;
#endif
#endif
}