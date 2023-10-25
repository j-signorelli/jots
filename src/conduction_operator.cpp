#include "conduction_operator.hpp"

using namespace std;


/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = m^{-1}(-Ku + boundary_terms)
 *
 *  where u is the vector representing the temperature, m is the mass matrix,
 *  and K is the stiffness matrix
 *
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
ConductionOperator::ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, const MaterialProperty* rho_prop, const MaterialProperty* C_prop, const MaterialProperty* k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref, const double& tf_ref)
:  TimeDependentOperator(f.GetTrueVSize(), t_ref),
   tf(tf_ref),
   time(t_ref),
   dt(dt_ref), 
   rho_C(rho_prop->GetLocalValue(0), C_prop->GetCoeffRef()),
   fespace(f),
   ode_solver(nullptr),
   expl_solver(nullptr),
   impl_solver(nullptr),
   m(nullptr), 
   k(nullptr),
   b(nullptr),
   M_full(nullptr),
   K_full(nullptr),
   b_vec(height),
   M(nullptr),
   K(nullptr),
   A(nullptr),
   M_e(nullptr),
   K_e(nullptr),
   A_e(nullptr),
   rhs(height),
   mass_updated(false)
{

   PreprocessSolver(in_config);

   PreprocessBCs(in_config, in_bcs, all_bdr_attr_markers);
   
   PreprocessMass();

   PreprocessStiffness(k_prop);

   //----------------------------------------------------------------
   // Instantiate ODESolver
   ode_solver = Factory::GetODESolver(in_config.GetTimeSchemeLabel());
   
   // Initialize ODESolver with THIS operator
   ode_solver->Init(*this);
}

void ConductionOperator::PreprocessSolver(const Config& in_config)
{  
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


}

void ConductionOperator::PreprocessBCs(const Config& in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers)
{

   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(in_config.GetBCCount());
   dbc_bdr = 0; // Start w/ all attributes set to non-essential = 0


   // Create linear form for Neumann BCs
   b = new ParLinearForm(&fespace);

   for (int i = 0; i < in_config.GetBCCount(); i++)
   {
      // Add Neumann BCs to linear form b
      if (!in_bcs[i]->IsEssential())
      {
         b->AddBoundaryIntegrator(new BoundaryLFIntegrator(in_bcs[i]->GetCoeffRef()), all_bdr_attr_markers[i]);// Add boundary integrator to boundary
      }
      else
      {
         // Update DBC list if essential
         dbc_bdr[i] = 1;
      }
   }

   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // Initialize Neumann linear form
   UpdateNeumann();
}

void ConductionOperator::PreprocessMass()
{   
   // Prepare mass matrix
   // Assemble parallel bilinear form for mass matrix
   m = new ParBilinearForm(&fespace);

   // Add the domain integral with included rho*Cp coefficient everywhere
   m->AddDomainIntegrator(new MassIntegrator(rho_C));

   // Initialize mass matrix data structures
   ReassembleMass();
}

void ConductionOperator::PreprocessStiffness(const MaterialProperty* k_prop)
{  

   // Assemble the parallel bilinear form for stiffness matrix
   k = new ParBilinearForm(&fespace);

   // Add domain integrator to the bilinear form with the cond_model coeff
   k->AddDomainIntegrator(new DiffusionIntegrator(k_prop->GetCoeffRef()));
   
   // Initialize stiffness data structures
   ReassembleStiffness();
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

void ConductionOperator::UpdateNeumann()
{
   b->Assemble();
   b->ParallelAssemble(b_vec);
}

void ConductionOperator::CalculateRHS(const Vector &u) const
{   

   // Complete multiplication of HypreParMatrix with primal vector for dual vector z
   K_full->Mult(u, rhs);

   rhs.Neg(); // z = -z

   // Add Neumann vector term
   rhs.Add(1, b_vec);

   // Now have RHS without "eliminated" essential DOFs
   // ^this must be done depending on LHS
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
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

bool ConductionOperator::IsNotComplete() const
{
   return !(time > tf || abs(time-tf) < TIME_TOLERANCE);
}

void ConductionOperator::Iterate(Vector& u)
{
   // Step in time
   ode_solver->Step(u, time, dt);
}

void ConductionOperator::ProcessMatPropUpdate(MATERIAL_PROPERTY mp)
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

ConductionOperator::~ConductionOperator()
{
   delete ode_solver;
   delete expl_solver;
   delete impl_solver;

   delete m;
   delete k;
   delete b;

   delete M_full;
   delete K_full;
   
   delete M;
   delete K;
   delete A;

   delete M_e;
   delete K_e;
   delete A_e;
}
