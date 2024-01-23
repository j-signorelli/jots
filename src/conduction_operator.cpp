#include "conduction_operator.hpp"

using namespace std;


double ConductionOperator::AOverBCCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return A->Eval(T, ip)/(B->Eval(T, ip)*C->Eval(T,ip));
}

double ConductionOperator::dAOverBCCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return dA->Eval(T, ip)/(B->Eval(T, ip)*C->Eval(T,ip)) 
           - A->Eval(T, ip)*dB->Eval(T,ip)/(pow(B->Eval(T, ip),2.0)*C->Eval(T,ip));
           - A->Eval(T, ip)*dC->Eval(T,ip)/(pow(C->Eval(T, ip),2.0)*B->Eval(T,ip))
}

ReducedSystemOperatorA::ReducedSystemOperatorA(Operator* K_, ParNonlinearForm* B_, ParNonlinearForm* N_)
: Operator(K_->Height()),
  K(K_), 
  B(B_),
  N(N_),
  N_vec(nullptr),
  Jacobian(nullptr)
{

}

ReducedSystemOperatorA::ReducedSystemOperatorA(Operator* K_, const Vector& N_vec_)
: Operator(K_->Height()),
  K(K_),
  B(nullptr),
  N(nullptr),
  N_vec(&N_vec_),
  Jacobian(nullptr)
{

}

void ReducedSystemOperatorA::Mult(const Vector &u, Vector &y) const
{
    K.Mult(u,y);

    if (B)
        B.AddMult(u,y)

    if (N)
        N.AddMult(u,y)
    else
        y += b_vec;// Else linear Neumann term -- just add
}

Operator& ReducedSystemOperatorA::GetGradient(const Vector& u) const
{
    delete Jacobian;
    
}


ConductionOperator::ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref)
:  TimeDependentOperator(f.GetTrueVSize(), t_ref),
   JOTSIterator(f, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()),
   time(t_ref),
   dt(dt_ref), 
   ode_solver(nullptr),
   expl_solver(nullptr),
   impl_solver(nullptr),
   newton_solver(f.GetComm(), true),
   one(1),
   zero(0),
   diffusitivity(k.GetCoeffRef(), rho.GetCoeffRef(), C.GetCoeffRef()),
   g_over_rhoC(neumann_coeff, rho.GetCoeffRef(), C.GetCoeffRef()),
   d_diffusivity(k.GetCoeffRef(), k.GetDCoeffRef(), rho.GetCoeffRef(), rho.GetDCoeffRef(), C.GetCoeffRef(), C.GetDCoeffRef()),
   dg_over_rhoC(neumann_coeff, zero, rho.GetCoeffRef(), rho.GetDCoeffRef(), C.GetCoeffRef(), C.GetDCoeffRef()),
   M(&f),
   K(nullptr),
   B(nullptr),
   N(nullptr),
   A(nullptr),
   R(&f),
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

    // Prepare mass bilinear form
    M.AddDomainIntegrator(new MassIntegrator());
	M.Assemble(0); // keep sparsity pattern of all the same!
    M.Finalize(0);

    // Stiffness matrix:
    if (k.IsConstant() && rho.IsConstant() && C.IsConstant())
    {
        // If all material properties are constant, then setup K as HypreParMatrix
        ParBilinearForm K_temp(diffusivity);
        K_temp.AddDomainIntegrator(new DiffusionIntegrator(diffusivity));
        K_temp.Assemble(0);
        K_temp.Finalize(0);
        K = K_temp.ParallelAssemble();// Ensure that K = PtAP
    }
    else
    {
        // Otherwise, prepare as ParNonlinearForm
        K = new ParNonlinearForm(&f);
        ParNonlinearForm* K_temp = dynamic_cast<ParNonlinearForm*>(K);
        K_temp->AddDomainIntegrator(new JOTSNonlinearDiffusionIntegrator(&f, diffusivity, d_diffusivity));
        K_temp->SetEssentialTrueDofs(ess_tdof_list);
    }
	
    if (!rho.IsConstant() || !C.IsConstant())
    {
        // Convection-type term for NL problems
        // If gradients of rho or C may exist, need to include this term
        B = new ParNonlinearForm(&f);
        B->AddDomainIntegrator(new JOTSNonlinearConvectionIntegrator(&f, beta_grad_u, dbeta_grad_u));
        B->SetEssentialTrueDofs(ess_tdof_list);

        // Neumann term
        N = new ParNonlinearForm(&f);
        N->AddBoundaryIntegrator(new JOTSNonlinearNeumannIntegrator(neumann, d_neumann));
        N->SetEssentialTrueDofs(ess_tdof_list);

        // Prepare ReducedSystemOperatorA with nonlinear Neumann term
        A = new ReducedSystemOperatorA(K, B, N);
    }
    else
    {
        // else leave B as null
        // Re-initialize Neumann term with division of input heat flux by rho*C
        delete b;
        b  = new ParLinearForm(fespace);
        b.AddBoundaryIntegrator(new BoundaryLFIntegrator(neumann));
        UpdateNeumann();
        // Prepare ReducedSystemOperatorA with linear Neumann term
        A = new ReducedSystemOperatorA(K, b_vec);
    }

	//----------------------------------------------------------------
	double abs_tol = in_config.GetAbsTol();
	double rel_tol = in_config.GetRelTol();
	int max_iter = in_config.GetMaxIter();
	//----------------------------------------------------------------
	// Prepare explicit solver

	expl_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());

	// Set up the solver for Mult
	expl_solver->iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false

	expl_solver->SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
	expl_solver->SetAbsTol(abs_tol); // Sets "absolute tolerance" of iterative solver
	expl_solver->SetMaxIter(max_iter); // Sets maximum number of iterations
	expl_solver->SetPrintLevel(0); // Print all information about detected issues

	expl_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel())); // Set type of preconditioning (relaxation type) 
	expl_solver->SetPreconditioner(expl_prec); // Set preconditioner to matrix inversion solver
	
	//----------------------------------------------------------------
	// Prepare implicit solver
	impl_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());

	// Set up solver for ImplicitSolve
	impl_solver->iterative_mode = false;
	impl_solver->SetRelTol(rel_tol);
	impl_solver->SetAbsTol(abs_tol);
	impl_solver->SetMaxIter(max_iter);
	impl_solver->SetPrintLevel(0);
	impl_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
	impl_solver->SetPreconditioner(impl_prec);

	//----------------------------------------------------------------
	// Prepare the NewtonSolver (for implicit)

	//----------------------------------------------------------------
	// Instantiate ODESolver
	ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());

	// Initialize ODESolver with THIS operator
	ode_solver->Init(*this);
}

void ConductionOperator::ReassembleMass()
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

void ConductionOperator::ReassembleStiffness()
{    

   delete K_full;
   delete K;
   delete K_e;

   K_full = nullptr;
   K = nullptr;
   K_e = nullptr;


   k->Update(); // delete old data (M and M_e) if conductivity changed

   k->Assemble(0);
   k->Finalize(0);

   // Create full stiffness matrix w/o removed essential DOFs
   K_full = k->ParallelAssemble();
   K = new HypreParMatrix(*K_full);

   // Create stiffness matrix w/ removed essential DOFs
   K_e = k->ParallelEliminateTDofs(ess_tdof_list, *K);


}

void ConductionOperator::CalculateRHS(const Vector &u) const
{   

   // Complete multiplication of HypreParMatrix with primal vector for dual vector z
   K_full->Mult(u, rhs);

   rhs.Neg(); // z = -z

   // Add Neumann vector term
   rhs.Add(1, b_vec_full);

   // Now have RHS without "eliminated" essential DOFs
   // ^this must be done depending on LHS
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
	// Compute:
	//    du/dt = M^{-1}(-Ku + b)
	// for du_dt

	// Evaluate action of nonlinear form for given u
	k.Mult(u, rhs);

	rhs.Neg(); // z = -z

	// Add Neumann term
	rhs.Add(1, b_vec);

    // Recompute M
    
    delete M;
    m->FormSystemMatrix

    // Reset explicit operator
    expl_solver->SetOperator(*M);
 
	// Calculate RHS pre-essential update
	CalculateRHS(u);
	
	du_dt.SetSubVector(ess_tdof_list, 0.0); // **Assume essential DOFs du_dt = 0.0

	// Apply elimination of essential BCs to rhs vector
	EliminateBC(*M, *M_e, ess_tdof_list, du_dt, rhs);

	// Solve for M^-1 * rhs
	expl_solver->Mult(rhs, du_dt);   

}


void ConductionOperator::ImplicitSolve(const double dt,
									   const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep

   // Deallocate any previous A, A_e
   delete A;
   delete A_e;

   // Calculate RHS pre-essential update
   CalculateRHS(u);

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

void ConductionOperator::Iterate(Vector& u)
{
   // Step in time
   ode_solver->Step(u, time, dt);
}

ConductionOperator::~ConductionOperator()
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
