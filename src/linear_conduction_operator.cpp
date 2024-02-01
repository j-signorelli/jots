#include "conduction_operator.hpp"

using namespace std;


/** After spatial discretization, the linear conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku + N)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the stiffness matrix
 *
 *  Class LinearConductionOperator represents the right-hand side of the above ODE.
 */
LinearConductionOperator::LinearConductionOperator(const Config &in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty &rho_prop, const MaterialProperty &C_prop, const MaterialProperty &k_prop, ParFiniteElementSpace &f, double &t_ref, double &dt_ref)
:  TimeDependentOperator(f.GetTrueVSize(), t_ref),
   JOTSIterator(f, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()),
   time(t_ref),
   dt(dt_ref), 
   rho_C(rho_prop.GetCoeffRef(), C_prop.GetCoeffRef()),
   ode_solver(nullptr),
   expl_solver(nullptr),
   impl_solver(nullptr),
   M(&f), 
   K(&f),
   M_full(nullptr),
   K_full(nullptr),
   M(nullptr),
   K(nullptr),
   A(nullptr),
   M_e(nullptr),
   K_e(nullptr),
   A_e(nullptr),
   rhs(height),
   mass_updated(false)
{
    double abs_tol = in_config.GetAbsTol();
    double rel_tol = in_config.GetRelTol();
    int max_iter = in_config.GetMaxIter();


    // Prepare explicit solver
    expl_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    expl_solver->iterative_mode = false;
    expl_solver->SetRelTol(rel_tol);
    expl_solver->SetAbsTol(abs_tol);
    expl_solver->SetMaxIter(max_iter);
    expl_solver->SetPrintLevel(0);
    expl_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel())); 
    expl_solver->SetPreconditioner(expl_prec);

    // Prepare implicit solver
    impl_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    impl_solver->iterative_mode = false;
    impl_solver->SetRelTol(rel_tol);
    impl_solver->SetAbsTol(abs_tol);
    impl_solver->SetMaxIter(max_iter);
    impl_solver->SetPrintLevel(0);
    impl_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
    impl_solver->SetPreconditioner(impl_prec);
    
    // Initialize mass
    M.AddDomainIntegrator(new MassIntegrator(rho_C));
    ReassembleMass();

    // Initialize stiffness
    K.AddDomainIntegrator(new DiffusionIntegrator(k_prop->GetCoeffRef()));
    ReassembleStiffness();

    // Instantiate ODESolver
    ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());
    
    // Initialize ODESolver with THIS operator
    ode_solver->Init(*this);
}

void LinearConductionOperator::ReassembleM()
{
    delete M_full;
    delete M;
    delete M_e;

    M_full = nullptr;
    M = nullptr;
    M_e = nullptr;

    m->Update(); // delete old data (M and M_e) if anything changed

    m->Assemble(0); // keep sparsity pattern of m and k the same
    m->Finalize(0);

    // Create full mass matrix w/o removed essential DOFs
    M_full = m->ParallelAssemble();
    M = new HypreParMatrix(*M_full);
    
    // Create mass matrix w/ removed essential DOFs + save eliminated portion
    M_e = m->ParallelEliminateTDofs(ess_tdof_list, *M);

    // Flag mass updated (called if Mult called)
    mass_updated = true;
}

void LinearConductionOperator::ReassembleK()
{    

    delete K_full;
    delete K;
    delete K_e;

    K_full = nullptr;
    K = nullptr;
    K_e = nullptr;


    k.Update(); // delete old data (M and M_e) if conductivity changed

    k.Assemble(0);
    k.Finalize(0);

    // Create full stiffness matrix w/o removed essential DOFs
    K_full = K.ParallelAssemble();
    K = new HypreParMatrix(*K_full);

    // Create stiffness matrix w/ removed essential DOFs
    K_e = K.ParallelEliminateTDofs(ess_tdof_list, *K);


}

void LinearConductionOperator::ReassembleA()
{

}

void LinearConductionOperator::CalculateRHS(const Vector &u, Vector &y) const
{   

    // Complete multiplication of HypreParMatrix with primal vector for dual vector z
    K_full->Mult(u, y);

    y.Neg(); // z = -z

    // Add Neumann vector term
    y.Add(1, b_vec_full);
}

void LinearConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du/dt = M^{-1}(-Ku + boundary_terms)
    // for du_dt

    // If mass matrix updated, reset operator
    if (mass_updated)
    {
        expl_solver->SetOperator(*M);
        mass_updated = false;
    }
    
    // Calculate RHS pre-essential update
    CalculateRHS(u, rhs);

    du_dt.SetSubVector(ess_tdof_list, 0.0); // **Assume essential DOFs du_dt = 0.0

    // Apply elimination of essential BCs to rhs vector
    EliminateBC(*M, *M_e, ess_tdof_list, du_dt, rhs);

    // Solve for M^-1 * rhs
    expl_solver->Mult(rhs, du_dt);   

}


void LinearConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
    // Solve the equation:
    //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
    // for du_dt, where K is linearized by using u from the previous timestep

    // Deallocate any previous A, A_e
    delete A;
    delete A_e;

    // Calculate RHS pre-essential update
    CalculateRHS(u, rhs);

    // Calculate LHS w/o essential DOFs eliminated at first
    A = Add(1.0, *M, dt, *K);

    impl_solver->SetOperator(*A);

    // Calculate LHS eliminated part using already calculated M_e's from the individual bilinear forms
    A_e = Add(1.0, *M_e, dt, *K_e);

    du_dt.SetSubVector(ess_tdof_list, 0.0); // **Assume essential DOFs du_dt = 0.0
                                            // Could this above be the issue^^?

    // Eliminate BCs from RHS given LHS
    A->EliminateBC(*A_e, ess_tdof_list, du_dt, rhs);

    // Solve for du_dt
    impl_solver->Mult(rhs, du_dt);
    }

void LinearConductionOperator::Iterate(Vector& u)
{
    // Step in time
    ode_solver->Step(u, time, dt);
}

void LinearConductionOperator::ProcessMatPropUpdate(MATERIAL_PROPERTY mp)
{
    switch (mp)
    {
        case MATERIAL_PROPERTY::DENSITY:
        case MATERIAL_PROPERTY::SPECIFIC_HEAT:
            ReassembleMass();
            break;
        case MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY:
            ReassembleStiffness();
            break;
    }
    }

LinearConductionOperator::~LinearConductionOperator()
{
    delete ode_solver;
    delete expl_solver;
    delete impl_solver;

    delete m;
    delete k;

    delete M_full;
    delete K_full;
    
    delete M;
    delete K;
    delete A;

    delete M_e;
    delete K_e;
    delete A_e;
}