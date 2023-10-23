#include "simulations.hpp"

using namespace mfem;
using namespace std;

const double UnsteadyHeatSimulation::TIME_TOLERANCE = 1e-14;

Simulation::Simulation(const std::string in_name, const mfem::ParGridFunction* u_0)
: solution_name(in_name)
{
    u_0->GetTrueDofs(u);
}

UnsteadyHeatSimulation::UnsteadyHeatSimulation(const mfem::ParGridFunction* u_0, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty* const* mat_props, ParFiniteElementSpace &f, double& in_time, double& in_dt)
: Simulation("Temperature", u_0),
  time(in_time),
  dt(in_dt),
  tf(in_config.GetFinalTime())
{
    // Instantiate ODESolver
    ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());

    // Instantiate ConductionOperator, sending all necessary parameters
    oper = new ConductionOperator(in_config, in_bcs, all_bdr_attr_markers, mat_props[MATERIAL_PROPERTY::DENSITY], mat_props[MATERIAL_PROPERTY::SPECIFIC_HEAT], mat_props[MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY], *f, time);
    
    // Initialize ODESolver with operator
    ode_solver->Init(*oper);

}

bool UnsteadyHeatSimulation::IsRunning()
{
    return !(time > tf || abs(time-tf) < TIME_TOLERANCE);
}

void UnsteadyHeatSimulation::Iterate()
{
    ode_solver->Step(u, time, dt);
}

UnsteadyHeatSimulation::~UnsteadyHeatSimulation()
{
    delete ode_solver;
    delete oper;
}