#include "steady_conduction_operator.hpp"

SteadyConductionOperator::SteadyConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& k_prop, ParFiniteElementSpace& f_)
: JOTSIterator(f_, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()),
  k(&f_),
  lin_solver(nullptr),
  newton(f_.GetComm())
{
    // Set linear solver + preconditioner
    lin_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    lin_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
    lin_solver->SetPreconditioner(lin_prec);
    lin_solver->SetAbsTol(in_config.GetAbsTol());
    lin_solver->SetRelTol(in_config.GetRelTol());
    lin_solver->SetMaxIter(in_config.GetMaxIter());

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

}

void SteadyConductionOperator::Iterate(Vector& u)
{
    // Solve
    newton.Mult(b, u);
}