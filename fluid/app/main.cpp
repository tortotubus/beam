#include <elff/fluid/Options.H>
#include <elff/fluid/Solver.H>

#include <AMReX.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ELFF::Fluid::Solver solver(ELFF::Fluid::SolverOptions::read());
        solver.initData();
        solver.evolve();
    }
    amrex::Finalize();
    return 0;
}
