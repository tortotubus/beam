

#include "embed.h"
#include "grid/octree.h"

#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
// #include "view.h"

face vector muv[];

int maxlevel = 10;
double D = 1. [1];
double L = 1. [1];
double U0 = 1.;
double Re = 300.;
int cylinder_endcaps = 1;
int cylinder_axis = 2; // 0: x-axis, 1: y-axis, 2: z-axis
// double theta = 0.;
// double phi = 0.;

char * base_path = "static_cylinder_3d_embed_output";

#define R (D / 2.)
#define CYLINDER_RADIAL_PHI                                                \
  (cylinder_axis == 0 ? sqrt(sq(y) + sq(z)) - R :                          \
   cylinder_axis == 1 ? sqrt(sq(x) + sq(z)) - R :                          \
                        sqrt(sq(x) + sq(y)) - R)
#define CYLINDER_AXIAL_PHI                                                 \
  (cylinder_axis == 0 ? fabs(x) - L / 2. :                                 \
   cylinder_axis == 1 ? fabs(y) - L / 2. :                                 \
                        fabs(z) - L / 2.)
#define CYLINDER_PHI                                                       \
  (cylinder_endcaps ? union(CYLINDER_RADIAL_PHI, CYLINDER_AXIAL_PHI)       \
                    : CYLINDER_RADIAL_PHI)

int main() {
  init_grid(64);
  size(16. * D);
  origin(-3. * D, -L0 / 2., -L0 / 2.);
  mu = muv;
  run();
}

event properties(i++) { foreach_face() muv.x[] = fm.x[] * D * U0 / Re; }

u.n[left] = dirichlet(U0);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

event init(t = 0) {
  if (cylinder_endcaps) {
    refine(CYLINDER_RADIAL_PHI < 0.2 * R && CYLINDER_AXIAL_PHI < 0.1 * L &&
           level < maxlevel);
    solid(cs, fs, CYLINDER_PHI);
  } else {
    refine(CYLINDER_RADIAL_PHI < 0.2 * R && level < maxlevel);
    solid(cs, fs, CYLINDER_PHI);
  }
  foreach () {
    u.x[] = cs[] ? U0 : 0.;
  }
}

event logfile(i++) fprintf(stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event geometry(t = 0) {
  if (pid() == 0)
    assert(!create_path(base_path));

  char fname[256];
  snprintf(fname, sizeof(fname), "%s/cylinder_facets.dat", base_path);
  FILE *fp = fopen(fname, "w");
  if (fp) {
    output_facets(cs, fp, fs);
    fclose(fp);
  }
}

event movies(t += 0.25; t <= 60) {
  scalar cylinder_phi[], l2[], vyz[];
  foreach()
    cylinder_phi[] = CYLINDER_PHI;
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] -
    u.z[0,-1]))/(2.*Delta);
  lambda2 (u, l2);

//   view (fov = 11.44, quat = {0.072072,0.245086,0.303106,0.918076},
// 	tx = -0.307321, ty = 0.22653, bg = {1,1,1},
// 	width = 802, height = 634);
//   draw_vof ("cs", "fs");
//   isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
// 	      linear = true, map = cool_warm);
//   save ("movie.mp4");
#if TREE
  output_hdf_htg({cs, cylinder_phi, l2, vyz, p}, {u}, base_path);
#else
  output_hdf_imagedata({cs, cylinder_phi, l2, vyz, p}, {u}, base_path);
#endif
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt(i++) {
  astats s =
      adapt_wavelet({cs, u}, (double[]){1e-2, 0.02, 0.02, 0.02}, maxlevel, 4);
  fprintf(stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
