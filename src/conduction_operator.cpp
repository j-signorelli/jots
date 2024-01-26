#include "conduction_operator.hpp"

using namespace std;

ReducedSystemOperatorA::ReducedSystemOperatorA(const Operator* K_, const ParNonlinearForm* B_, const ParNonlinearForm* N_)
: Operator(K_->Height()),
  K_mat(dynamic_cast<const HypreParMatrix*>(K_)),
  K(dynamic_cast<const ParNonlinearForm*>(K_)), 
  B(B_),
  N(N_),
  N_vec(nullptr),
  Jacobian(nullptr)
{

}

ReducedSystemOperatorA::ReducedSystemOperatorA(const Operator* K_, const Vector* N_vec_)
: Operator(K_->Height()),
  K_mat(dynamic_cast<const HypreParMatrix*>(K_)),
  K(dynamic_cast<const ParNonlinearForm*>(K_)),
  B(nullptr),
  N(nullptr),
  N_vec(N_vec_),
  Jacobian(nullptr)
{

}

void ReducedSystemOperatorA::Mult(const Vector &u, Vector &y) const
{
    if (K)
        K->Mult(u,y);
    else
        K_mat->Mult(u,y);

    y.Neg();

    if (B)
    {
        B->AddMult(u,y);
        y.Neg();
        N->AddMult(u,y);
    }
    else
    {
        y += *N_vec;// Else just linear Neumann term -- just add
    }
}

Operator& ReducedSystemOperatorA::GetGradient(const Vector& u) const
{
    delete Jacobian;
    Jacobian = nullptr;
    // If K is linear (ie: k, rho, C all constant), then Jacobian is zero
    // If above true, then operator is a HypreParMatrix
    if (K)
    {
        // Default Operator type for ParNonlinearForm::GetGradient are &HypreParMatrix
        // So can static_cast all safely
        HypreParMatrix &K_grad = static_cast<HypreParMatrix&>(K->GetGradient(u));
        Jacobian = new HypreParMatrix(K_grad); // Create copy
        *Jacobian *= -1;
        if (B) // if rho(u) or C(u)
        {
            Jacobian->Add(-1, static_cast<HypreParMatrix&>(B->GetGradient(u)));
            Jacobian->Add(1, static_cast<HypreParMatrix&>(N->GetGradient(u)));
        }
    }
    else
    {
        Jacobian = new HypreParMatrix(*K_mat); // Just create copy of K_mat
        *Jacobian *= -1;
    }

    return *Jacobian;

}


void ReducedSystemOperatorR::Mult(const Vector &k, Vector &y) const
{
    add(*u_n, *dt, k, z);
    A.Mult(z, y);
    y.Neg();

    M_mat.AddMult(k, y);
    y.SetSubVector(ess_tdof_list, 0.0); // Pertinent if linear
}

Operator& ReducedSystemOperatorR::GetGradient(const Vector &k) const
{
    delete Jacobian;
    Jacobian = nullptr;
    add(*u_n, *dt, k, z);
    Jacobian = Add(1.0, M_mat, -1.0*(*dt), static_cast<HypreParMatrix&>(A.GetGradient(z)));
    HypreParMatrix* Je = Jacobian->EliminateRowsCols(ess_tdof_list);
    delete Je; // ^Not sure if above is necessary -- TODO
    return *Jacobian;
}

double ConductionOperator::DiffusivityCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return k.Eval(T, ip)/(rho.Eval(T, ip)*C.Eval(T,ip));
}

double ConductionOperator::dDiffusivityCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    double k_val = k.Eval(T, ip);
    double dk_val = dk.Eval(T, ip);
    double rho_val = rho.Eval(T, ip);
    double drho_val = drho.Eval(T, ip);
    double C_val = C.Eval(T, ip);
    double dC_val = dC.Eval(T, ip);
    
    return dk_val/(rho_val*C_val) 
           - k_val*drho_val/(rho_val*rho_val*C_val)
           - k_val*dC_val/(C_val*C_val*rho_val);
}

double ConductionOperator::NeumannCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    return g.Eval(T, ip)/(rho.Eval(T, ip)*C.Eval(T,ip));
}

double ConductionOperator::dNeumannCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    double g_val = g.Eval(T, ip);

    double rho_val = rho.Eval(T, ip);
    double drho_val = drho.Eval(T, ip);
    double C_val = C.Eval(T, ip);
    double dC_val = dC.Eval(T, ip);
    
    return - g_val*drho_val/(rho_val*rho_val*C_val)
           - g_val*dC_val/(C_val*C_val*rho_val);
}

double ConductionOperator::BetaCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    double k_val = k.Eval(T, ip);

    double rho_val = rho.Eval(T, ip);
    double drho_val = drho.Eval(T, ip);
    double C_val = C.Eval(T, ip);
    double dC_val = dC.Eval(T, ip);
    return k_val*(-drho_val/(rho_val*rho_val*C_val)
                  -dC_val/(C_val*C_val*rho_val));
}

double ConductionOperator::dBetaCoefficient::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
    double k_val = k.Eval(T, ip);
    double dk_val = dk.Eval(T, ip);
    double rho_val = rho.Eval(T, ip);
    double drho_val = drho.Eval(T, ip);
    double d2rho_val = d2rho.Eval(T, ip);
    double C_val = C.Eval(T, ip);
    double dC_val = dC.Eval(T, ip);
    double d2C_val = d2C.Eval(T, ip);

    return dk_val*(-drho_val/(rho_val*rho_val*C_val)
                  -dC_val/(C_val*C_val*rho_val))
            + k_val*(-d2C_val/(C_val*C_val*rho_val)
                    +  2*dC_val*drho_val/(C_val*C_val*rho_val*rho_val)
                    + 2*dC_val*dC_val/(C_val*C_val*C_val*rho_val)
                    - d2rho_val/(C_val*rho_val*rho_val)
                    + 2*drho_val*drho_val/(C_val*rho_val*rho_val*rho_val));

}


ConductionOperator::ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double& time_, double& dt_)
:  TimeDependentOperator(f.GetTrueVSize(), time_),
   JOTSIterator(f, in_bcs, all_bdr_attr_markers, in_config.GetBCCount()),
   time(time_),
   dt(dt_),
   ode_solver(nullptr),
   lin_solver(nullptr),
   newton(f.GetComm()),
   diffusivity(k_prop, rho_prop, C_prop),
   d_diffusivity(k_prop, rho_prop, C_prop),
   g_over_rhoC(neumann_coeff, rho_prop, C_prop),
   dg_over_rhoC(neumann_coeff, rho_prop, C_prop),
   beta(k_prop, rho_prop, C_prop),
   d_beta(k_prop, rho_prop, C_prop),
   M(&f),
   K(nullptr),
   B(&f),
   N(&f),
   A(nullptr),
   R(nullptr),
   rhs(height)
{

    // Prepare mass bilinear form
    M.AddDomainIntegrator(new MassIntegrator());
	M.Assemble(0); // keep sparsity pattern of all the same!
    M.Finalize(0);
    M.FormSystemMatrix(ess_tdof_list, M_mat);

    // Stiffness matrix:
    if (k_prop.IsConstant() && rho_prop.IsConstant() && C_prop.IsConstant())
    {
        // If all material properties are constant, then setup K as HypreParMatrix
        // ^Avoids reassembly
        ParBilinearForm K_temp(&f);
        K_temp.AddDomainIntegrator(new DiffusionIntegrator(diffusivity));
        K_temp.Assemble(0);
        K_temp.Finalize(0);
        K = K_temp.ParallelAssemble();// Ensure that K = PtAP
    }
    else
    {
        // Otherwise, prepare as K ParNonlinearForm
        K = new ParNonlinearForm(&f);
        ParNonlinearForm* K_temp = dynamic_cast<ParNonlinearForm*>(K);
        K_temp->AddDomainIntegrator(new JOTSNonlinearDiffusionIntegrator(&f, diffusivity, d_diffusivity));
        K_temp->SetEssentialTrueDofs(ess_tdof_list);

    }

    if (!rho_prop.IsConstant() || !C_prop.IsConstant())
    {
        // If NL Neumann term...
        // Convection-type term for NL problems
        // If gradients of rho or C may exist, need to include this term
        B.AddDomainIntegrator(new JOTSNonlinearConvectionIntegrator(&f, beta, d_beta));
        B.SetEssentialTrueDofs(ess_tdof_list);

        // Neumann term
        N.AddBoundaryIntegrator(new JOTSNonlinearNeumannIntegrator(g_over_rhoC, dg_over_rhoC));
        N.SetEssentialTrueDofs(ess_tdof_list);

        // Prepare ReducedSystemOperatorA with nonlinear Neumann term
        A = new ReducedSystemOperatorA(K, &B, &N);
    }
    else
    {   
        // Otherwise:
        // Re-initialize Neumann term with division of input heat flux by rho*C
        // ^Avoids any reassembly
        b.GetBLFI()->DeleteAll(); // Delete all current boundary linear form integrators
        b.AddBoundaryIntegrator(new BoundaryLFIntegrator(g_over_rhoC));
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
	// Prepare the NewtonSolver
    // Use same linear solver obj and settings
    R = new ReducedSystemOperatorR(M_mat, *A, ess_tdof_list);
    newton.iterative_mode = false;
    newton.SetOperator(*R);
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
	A->Mult(u, rhs);
    
    // Zero out essential tdof rows on RHS
	rhs.SetSubVector(ess_tdof_list, 0.0); // Pertinent if linear

	// Solve for M^-1 * rhs
    lin_solver->SetOperator(M_mat);
	lin_solver->Mult(rhs, k);   

}


void ConductionOperator::ImplicitSolve(const double dt,
									   const Vector &u, Vector &k)
{

    // Solve the equation:
    //    R(k_{n+1}) = Mk_{n+1} - A(u_n + k*dt, t_n) = 0
    // for k

    // Update R(k) w/ most recent state
    R->SetParameters(&u, &dt);
    rhs = 0.0;
    newton.Mult(rhs, k);
} 

void ConductionOperator::Iterate(Vector& u)
{
   // Step in time
   double dt_ = dt; // make copy of dt
   ode_solver->Step(u, time, dt);
   // Verify no change in dt
   MFEM_VERIFY(dt_ == dt, "JOTS does not support time integration shemes that implicitly change dt!");
   t = time;
}

ConductionOperator::~ConductionOperator()
{
   delete ode_solver;
   delete lin_solver;

   delete K;
   delete A;
   delete R;
}
