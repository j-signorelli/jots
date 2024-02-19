#include "linear_conduction_operator.hpp"

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
   M_mat(nullptr),
   K_mat(nullptr),
   A_mat(nullptr),
   rhs(height)
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
    expl_solver->SetPrintLevel(Factory::CreatePrintLevel(in_config.GetLinSolPrintLevel()));

    // Prepare implicit solver
    impl_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    impl_solver->iterative_mode = false;
    impl_solver->SetRelTol(rel_tol);
    impl_solver->SetAbsTol(abs_tol);
    impl_solver->SetMaxIter(max_iter);
    impl_solver->SetPrintLevel(0);
    impl_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
    impl_solver->SetPreconditioner(impl_prec);
    impl_solver->SetPrintLevel(Factory::CreatePrintLevel(in_config.GetLinSolPrintLevel()));
    
    // Initialize mass
    M.AddDomainIntegrator(new MassIntegrator(rho_C));
    ReassembleM();

    // Initialize stiffness
    K.AddDomainIntegrator(new DiffusionIntegrator(k_prop.GetCoeffRef()));
    ReassembleK();

    // Initialize A
    ReassembleA();

    // Instantiate ODESolver
    ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());
    
    // Initialize ODESolver with THIS operator
    ode_solver->Init(*this);
}

void LinearConductionOperator::ReassembleM()
{
    delete M_mat;
    M_mat = nullptr;

    M.Update(); // delete old data (M and M_e) if anything changed

    M.Assemble(0); // keep sparsity pattern of m and k the same
    M.Finalize(0);

    M_mat = M.ParallelAssemble();

    M_mat->EliminateBC(ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);

    // Reset operator for explicit solver
    expl_solver->SetOperator(*M_mat);
}

void LinearConductionOperator::ReassembleK()
{    
    delete K_mat;
    K_mat = nullptr;

    K.Update(); // delete old data (M and M_e) if conductivity changed

    K.Assemble(0);
    K.Finalize(0);
    

    K_mat = K.ParallelAssemble();
    // Notably, do NOT apply any elimination of BCs to K
}

void LinearConductionOperator::ReassembleA()
{
    delete A_mat;
    A_mat = nullptr;

     // Mixing of eliminated + non-eliminated BCs OK here as neither M nor K not multiplied by anything
    A_mat = Add(1.0, *M_mat, dt, *K_mat);

    A_mat->EliminateBC(ess_tdof_list, Operator::DiagonalPolicy::DIAG_ONE);

    impl_solver->SetOperator(*A_mat);
}

void LinearConductionOperator::CalculateRHS(const Vector &u, Vector &y) const
{   

    // Complete multiplication of K with primal tvector for dual tvector z
    K_mat->Mult(u, y);

    y.Neg(); // -Ku

    // Add Neumann vector term
    y.Add(1.0, b_vec);
}

void LinearConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
    // Compute:
    //    du/dt = M^{-1}(-Ku + boundary_terms)
    // for du_dt
    
    // Calculate RHS
    CalculateRHS(u, rhs);

    // Zero out essential DOF rows -- as we want du_p/dt = 0
    rhs.SetSubVector(ess_tdof_list, 0.0);

    // Solve for M^-1 * rhs
    expl_solver->Mult(rhs, du_dt);

}


void LinearConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
    // Solve the equation:
    //    M*du_dt = [-K(u + dt*du_dt) + N]
    // (M + dt*K)*du_dt = -Ku + N
    // for du_dt, where K is linearized by using u from the previous timestep

    // Calculate RHS pre-essential update
    CalculateRHS(u, rhs);

    // Zero out essential DOF rows -- as we want du_p/dt = 0
    rhs.SetSubVector(ess_tdof_list, 0.0);

    // Solve for du_dt
    impl_solver->Mult(rhs, du_dt);

}

void LinearConductionOperator::Iterate(Vector& u)
{
    // Step in time
    double dt_ = dt;
    ode_solver->Step(u, time, dt);
    MFEM_VERIFY(dt_ == dt, "JOTS does not support time integration shemes that implicitly change dt!");
    t = time;
}

void LinearConductionOperator::ProcessMatPropUpdate(MATERIAL_PROPERTY mp)
{
    switch (mp)
    {
        case MATERIAL_PROPERTY::DENSITY:
        case MATERIAL_PROPERTY::SPECIFIC_HEAT:
            ReassembleM();
            ReassembleA();
            break;
        case MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY:
            ReassembleK();
            ReassembleA();
            break;
    }
    }

LinearConductionOperator::~LinearConductionOperator()
{
    delete ode_solver;
    delete expl_solver;
    delete impl_solver;

    delete M_mat;
    delete K_mat;
    delete A_mat;
}