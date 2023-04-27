#include "solver_state.hpp"

using namespace mfem;
using namespace std;

SolverState::SolverState(mfem::ParFiniteElementSpace* f)
: it_num(0),
  time(0.0),
  dt(0.0),
  tf(0.0),
  fespace(f)
{
    T_gf = new ParGridFunction(fespace);
}

SolverState::~SolverState()
{
    delete T_gf;
}