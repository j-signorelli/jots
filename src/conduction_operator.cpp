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
    y.Neg();

    if (B)
    {
        B.AddMult(u,y)
        y.Neg();
        N.AddMult(u,y)
    }
    else
        y += b_vec;// Else just linear Neumann term -- just add

}

Operator& ReducedSystemOperatorA::GetGradient(const Vector& u) const
{
    delete Jacobian;

    // If K is linear (ie: k, rho, C all constant), then Jacobian is zero
    // If above true, then operator is a HypreParMatrix
    if (K.GetType() == Operator::Hypre_ParCSR)
    {
        MFEM_WARNING("ReducedSystemOperatorA::GetGradient is being called despite all terms being linear. Newton iterations should be set to 1." )
        Jacobian = new HypreParMatrix(dynamic_cast<HypreParMatrix>(K));
        Jacobian = 0;
    }
    else
    {
        Jacobian = new HypreParMatrix(K.GetGradient(u));
        Jacobian *= -1;
        if (B) // if rho(u) or C(u)
        {
            Jacobian -= B.GetGradient(u);
            Jacobian += N.GetGradient(u);
        }
    }

    return *Jacobian;

}


void ReducedSystemOperatorR::Mult(const Vector &k, Vector &y) const
{
    add(*u, dt, k, z);
    A.Mult(z, y);
    y.Neg();
    M.TrueAddMult(k, y); // Ensure y += (PtMP) k
    y.SetSubVector(ess_tdof_list, 0.0); // Pertinent if linear
}

Operator& ReducedSystemOperatorR::GetGradient(const Vector &u) const
{
    delete Jacobian;
    add(*u, dt, k, z);
    Jacobian = Add(1.0, ..., -dt, A.GetGradient(z));
    return *Jacobian;
}

ConductionOperator::ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref)
:  TimeDependentOperator(f.GetTrueVSize(), t_ref),
   JOTSIterator(f, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()),
   ode_solver(nullptr),
   lin_solver(nullptr),
   impl_solver(nullptr),
   newton_solver(f.GetComm(), true),
   diffusitivity(k.GetCoeffRef(), rho.GetCoeffRef(), C.GetCoeffRef()),
   g_over_rhoC(neumann_coeff, rho.GetCoeffRef(), C.GetCoeffRef()),
   d_diffusivity(k.GetCoeffRef(), k.GetDCoeffRef(), rho.GetCoeffRef(), rho.GetDCoeffRef(), C.GetCoeffRef(), C.GetDCoeffRef()),
   dg_over_rhoC(neumann_coeff, zero, rho.GetCoeffRef(), rho.GetDCoeffRef(), C.GetCoeffRef(), C.GetDCoeffRef()),
   M(&f),
   K(nullptr),
   B(nullptr),
   N(nullptr),
   A(nullptr),
   R(nullptr)
   rhs(height)
{

    // Prepare mass bilinear form
    M.AddDomainIntegrator(new MassIntegrator());
	M.Assemble(0); // keep sparsity pattern of all the same!
    M.Finalize(0);
    M.FormSystemMatrix(ess_tdof_list, Mmat);

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
        A = new ReducedSystemOperatorA(K, &b_vec);
    }

	//----------------------------------------------------------------
	double abs_tol = in_config.GetAbsTol();
	double rel_tol = in_config.GetRelTol();
	int max_iter = in_config.GetMaxIter();
	//----------------------------------------------------------------
	// Prepare linear solver
	lin_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());

	// Set up the linear solver
	lin_solver->iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false
	lin_solver->SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
	lin_solver->SetAbsTol(abs_tol); // Sets "absolute tolerance" of iterative solver
	lin_solver->SetMaxIter(max_iter); // Sets maximum number of iterations
	lin_solver->SetPrintLevel(0); // Print all information about detected issues
	lin_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel())); // Set type of preconditioning (relaxation type) 
	lin_solver->SetPreconditioner(lin_prec); // Set preconditioner to matrix inversion solver

	//----------------------------------------------------------------
	// Prepare the NewtonSolver (for implicit time integration)
    // Use same linear solver obj and settings
    R = new ReducedSystemOperatorR(M, A);
    newton.iterative_mode = false;
    newton.SetOperator(R);
    newton.AddMaterialProperty(k_prop); // Register k_prop to be updated with solution vector every Newton iteration
    newton.AddMaterialProperty(C_prop);
    newton.AddMaterialProperty(rho_prop);
    newton.SetSolver(*lin_solver);
    newton.SetPrintLevel(1);
    newton.SetAbsTol(in_config.GetNewtonAbsTol());
    newton.SetRelTol(in_config.GetNewtonRelTol());
    newton.SetMaxIter(in_config.GetNewtonMaxIter());

	//----------------------------------------------------------------
	// Instantiate ODESolver
	ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());

	// Initialize ODESolver with THIS operator
	ode_solver->Init(*this);
}


void ConductionOperator::Mult(const Vector &u, Vector &k) const
{
	// Compute:
	//    k_n = M^{-1}A(u_n,t_n)
	// for k

	// Evaluate action of A
	A.Mult(u, rhs);
    
    // Zero out essential tdof rows on RHS
	rhs.SetSubVector(ess_tdof_list, 0.0); // Pertinent if linear

	// Solve for M^-1 * rhs
    lin_solver->SetOperator(*M_mat);
	lin_solver->Mult(rhs, k);   

}


void ConductionOperator::ImplicitSolve(const double dt,
									   const Vector &u, Vector &k)
{
    // Solve the equation:
    //    R(k_{n+1}) = Mk_{n+1} - A(u_n + k*dt, t_n) = 0

    // Update R(k) w/ most recent state
    R.SetParameters(&u, &dt)
    Vector zero;
    newton.Mult(zero, du_dt);
} 

void ConductionOperator::Iterate(Vector& u)
{
   // Step in time
   ode_solver->Step(u, time, dt);
}

ConductionOperator::~ConductionOperator()
{
   delete ode_solver;
   delete lin_solver;

   delete K;
   delete B;
   delete N;
   delete A;
   delete R;
}
