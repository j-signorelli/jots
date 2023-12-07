#include "steady_conduction_operator.hpp"


using namespace mfem;

NonlinearJOTSDiffusionIntegrator::NonlinearJOTSDiffusionIntegrator(MaterialProperty& k_, ParFiniteElementSpace* fespace_)
: k(k_),
  fespace(*fespace_),
  u_gf(fespace_),
  u_coeff(&u_gf),
  dkdu_times_u(k.GetDCoeffRef(), u_coeff),
  diff(k.GetCoeffRef()),
  diff_d(dkdu_times_u)
{

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Update the coefficient k
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);

    // Very simply can just use DiffusionIntegrator::AssembleElementVector,
    // as the coefficient is re-evaluated every call
    diff.AssembleElementVector(el, Tr, elfun, elvect);

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    // Update both coefficients as both required to be updated
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);
    k.UpdateDCoeff(elfun, dofs);
    
    // Update u GridFunction associated with dkdu * u
    u_gf.SetSubVector(dofs, elfun);

    // Sum of usual diffusion term plus nonlinear term
    // Very simply can just use DiffusionIntegrator::AssembleElementMatrix,
    // as the coefficient is re-evaluated every call
    diff.AssembleElementMatrix(el, Tr, elmat);

    // Can use DiffusionIntegrator but with lambda = k'(u)u
    DenseMatrix elmat_d;
    diff_d.AssembleElementMatrix(el, Tr, elmat_d);

    elmat += elmat_d;


}

SteadyConductionOperator::SteadyConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& k_prop, ParFiniteElementSpace& f_)
: JOTSIterator(f_, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()), 
  k(&f_),
  lin_solver(nullptr)
{
    // Add nonlinear diffusion
    k.AddDomainIntegrator(new NonlinearJOTSDiffusionIntegrator(k_prop, &f_));
    k.SetEssentialTrueDofs(ess_tdof_list);

    
    double abs_tol = in_config.GetAbsTol();
    double rel_tol = in_config.GetRelTol();
    int max_iter = in_config.GetMaxIter();

    // Set linear solver + preconditioner
    lin_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    lin_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
    lin_solver->SetPreconditioner(lin_prec);

    // Set NewtonSolver
    newton.SetOperator(k);
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