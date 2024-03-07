#include "eq_linear_elastic_operator.hpp"

EquilibriumLinearElasticOperator::EquilibriumLinearElasticOperator(ParFiniteElementSpace& f_, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty &lambda_prop, MaterialProperty &mu_prop)
: JOTSIterator(f_, in_config, in_bcs, all_bdr_attr_markers, f_.GetParMesh()->bdr_attributes.Size()),
  a(&f_),
  A_mat(nullptr)
{
    lin_solver->iterative_mode = true;
    
    // Add elasticity integrator
    a.AddDomainIntegrator(new ElasticityIntegrator(lambda_prop.GetCoeffRef(), mu_prop.GetCoeffRef()));
    a.Assemble(0);
    a.Finalize(0);
    A_mat = a.ParallelAssemble();

    A_mat->EliminateBC(ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);
    lin_solver->SetOperator(*A_mat);

}

void EquilibriumLinearElasticOperator::Iterate(Vector& u)
{
    // Solve
    lin_solver->Mult(b_vec, u);
}
    