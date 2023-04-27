#include "solver_state.hpp"

using namespace mfem;
using namespace std;

SolverState::SolverState(mfem::ParGridFunction initial_gf, int in_it, double in_time, double in_dt, const double in_tf);
: it_num(in_int),
  time(in_time,
  dt(in_dt),
  tf(in_tf),
  fespace(*initial_gf.ParFESpace())
{
    initial_gf.SetTrueDofs(T);
}
