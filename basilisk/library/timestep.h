#define CKPT_TIMESTEP 1

struct CheckpointSidecar {
  int iter;
  int inext;
  double time;
  double dt;
  double tnext;
};

static struct CheckpointSidecar cpsc = {-1,-1,0.,0.,0.};

static int write_checkpoint_sidecar (const char* file,
                                     double time,
                                     double dt_,
                                     int iter_,
                                     int inext_,
                                     double tnext_) {
  FILE* fp = fopen (file, "wb");
  if (!fp)
    return -1;

  struct CheckpointSidecar s = {.iter = iter_,
                                .inext = inext_,
                                .time = time,
                                .dt = dt_,
                                .tnext = tnext_};

  int ok = fwrite (&s, sizeof (s), 1, fp) == 1 ? 0 : -1;
  fclose (fp);
  return ok;
}

static int read_checkpoint_sidecar (const char* file,
                                    struct CheckpointSidecar* s) {
  FILE* fp = fopen (file, "rb");
  if (!fp)
    return -1;

  int ok = fread (s, sizeof (*s), 1, fp) == 1 ? 0 : -1;
  fclose (fp);

  return ok;
}

double timestep (const face vector u, double dtmax)
{
  static double previous = 0.;

  if (cpsc.iter >= 0) {
    previous = cpsc.dt;
    dt = cpsc.dt;
    cpsc.iter = -1;
    return dt;
  }
  if (t == 0.) previous = 0.;
  dtmax /= CFL;

  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
      assert (fm.x[]);
      dt *= fm.x[];
      if (dt < dtmax) dtmax = dt;
    }

  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
