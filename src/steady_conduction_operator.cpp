#include "steady_conduction_operator.hpp"

SteadyConductionOperator::SteadyConductionOperator(ParFiniteElementSpace& f_, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& k_prop)
: JOTSIterator(f_, in_config, in_bcs, all_bdr_attr_markers, f_.GetParMesh()->bdr_attributes.Size()),
  k(&f_),
  newton(f_.GetComm())
{
    lin_solver->iterative_mode = true;

    // Add nonlinear diffusion
    k.AddDomainIntegrator(new JOTSNonlinearDiffusionIntegrator(&f_, k_prop.GetCoeffRef(), k_prop.GetDCoeffRef()));
    k.SetEssentialTrueDofs(ess_tdof_list);

    // Set NewtonSolver
    newton.SetOperator(k);
    newton.AddMaterialProperty(k_prop); // Register k_prop to be updated with solution vector every Newton iteration
    newton.SetSolver(*lin_solver);
    newton.SetPrintLevel(1);
    newton.SetAbsTol(in_config.GetNewtonAbsTol());
    newton.SetRelTol(in_config.GetNewtonRelTol());
    newton.SetMaxIter(in_config.GetNewtonMaxIter());
    newton.SetPrintLevel(Factory::CreatePrintLevel(in_config.GetNewtonPrintLevel()));

}

void SteadyConductionOperator::Iterate(Vector& u)
{
    // Solve
    newton.Mult(b_vec, u);
}