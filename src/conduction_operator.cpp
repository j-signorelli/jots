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
ConductionOperator::ConductionOperator(const Config* in_config, BoundaryCondition** in_bcs, const ConductivityModel* in_cond, ParFiniteElementSpace &f, double t_0)
   :  TimeDependentOperator(f.GetTrueVSize(), t_0),
      fespace(f),
      impl_solver(NULL),
      expl_solver(NULL),
      m(NULL), 
      k(NULL),
      b(NULL),
      M_full(NULL),
      K_full(NULL),
      b_vec(height),
      M(NULL),
      K(NULL),
      A(NULL),
      M_e(NULL),
      K_e(NULL),
      A_e(NULL),
      rhs(height),
      user_input(in_config),
      boundary_conditions(in_bcs),
      constant_cond_model(in_cond->IsConstant())
{

   PreprocessBCs();
   
   PreprocessStiffness(in_cond);

   PreprocessSolver();

}

void ConductionOperator::PreprocessBCs()
{

   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(user_input->GetBCCount());
   dbc_bdr = 0; // Start w/ all attributes set to non-essential = 0

   // Create bdr_attr_marker arrays and initialize dbc list
   all_bdr_attr_markers = new Array<int>[user_input->GetBCCount()];
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      Array<int> bdr_attr(user_input->GetBCCount());
      bdr_attr = 0;
      bdr_attr[i] = 1;
      all_bdr_attr_markers[i] = bdr_attr;

      // Update DBC list if essential
      if (boundary_conditions[i]->IsEssential())
         dbc_bdr[i] = 1;
   }

   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // Create linear form for Neumann BCs
   b = new ParLinearForm(&fespace);

   // Loop through all BCs
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {  
      // Initialize coefficients for all BCs
      boundary_conditions[i]->InitCoefficient();
   
      // Add Neumann BCs to linear form b
      if (! boundary_conditions[i]->IsEssential())
      {
         b->AddBoundaryIntegrator(new BoundaryLFIntegrator(*boundary_conditions[i]->GetCoeffPtr()), all_bdr_attr_markers[i]);// Add boundary integrator to boundary
      }
   }


   // Assemble b
   b->Assemble();

   // Set b_vec
   b->ParallelAssemble(b_vec);
}

void ConductionOperator::PreprocessStiffness(const ConductivityModel* in_cond)
{  

   // Assemble the parallel bilinear form for stiffness matrix
   k = new ParBilinearForm(&fespace);

   // Add domain integrator to the bilinear form with the cond_model coeff
   k->AddDomainIntegrator(new DiffusionIntegrator(in_cond->GetCoeffRef()));
   
   // Initialize stiffness data structures
   UpdateStiffness();
}

void ConductionOperator::PreprocessSolver()
{  
   double abs_tol = user_input->GetAbsTol();
   double rel_tol = user_input->GetRelTol();
   int max_iter = user_input->GetMaxIter();

   //----------------------------------------------------------------   
   // Prepare mass matrix
   // Assemble parallel bilinear form for mass matrix
   m = new ParBilinearForm(&fespace);

   // Add the domain integral with included rho*Cp coefficient everywhere
   ConstantCoefficient rhoCp(user_input->GetDensity()*user_input->GetCp());
   m->AddDomainIntegrator(new MassIntegrator(rhoCp));
   m->Assemble(0); // keep sparsity pattern of m and k the same
   m->Finalize(0);

   // Create full mass matrix w/o removed essential DOFs
   M_full = m->ParallelAssemble();
   M = new HypreParMatrix(*M_full);
   
   // Create mass matrix w/ removed essential DOFs + save eliminated portion
   M_e = m->ParallelEliminateTDofs(ess_tdof_list, *M);

   //----------------------------------------------------------------
   // Prepare explicit solver
   
   expl_solver = user_input->GetSolver(fespace.GetComm());

   // Set up the solver for Mult
   expl_solver->iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false

   expl_solver->SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
   expl_solver->SetAbsTol(abs_tol); // Sets "absolute tolerance" of iterative solver
   expl_solver->SetMaxIter(max_iter); // Sets maximum number of iterations
   expl_solver->SetPrintLevel(0); // Print all information about detected issues

   expl_prec.SetType(user_input->GetPrec()); // Set type of preconditioning (relaxation type) 
   expl_solver->SetPreconditioner(expl_prec); // Set preconditioner to matrix inversion solver

   expl_solver->SetOperator(*M); // Set operator M

   //----------------------------------------------------------------
   // Prepare implicit solver
   impl_solver = user_input->GetSolver(fespace.GetComm());
   
   // Prepare bilinear form for LHS
   //a = new ParBilinearForm(&fespace);

   // Set up solver for ImplicitSolve
   impl_solver->iterative_mode = false;
   impl_solver->SetRelTol(rel_tol);
   impl_solver->SetAbsTol(abs_tol);
   impl_solver->SetMaxIter(max_iter);
   impl_solver->SetPrintLevel(0);
   impl_prec.SetType(user_input->GetPrec());
   impl_solver->SetPreconditioner(impl_prec);
   

}

void ConductionOperator::PreprocessIteration(Vector &u)
{
   // Apply BCs
   ApplyBCs(u);

}

void ConductionOperator::ApplyBCs(Vector &u)
{

   // Update time-dependent coefficients with current time
   // Also project essential coefficients onto u
   ParGridFunction temp_u_gf(&fespace);
   temp_u_gf.SetFromTrueDofs(u);

   bool n_changed = false;
   bool d_changed = false;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      if (!boundary_conditions[i]->IsConstant()) // If not constant in time
      {
         // Update coefficients (could be preCICE calls, could be SetTime calls, etc.)
         boundary_conditions[i]->UpdateCoeff();
         n_changed = true;
      }

      if (boundary_conditions[i]->IsEssential())
      {
         // Project correct values on boundary for essential BCs
         d_changed = true;
         temp_u_gf.ProjectBdrCoefficient(*boundary_conditions[i]->GetCoeffPtr(), all_bdr_attr_markers[i]);
      }
   }

   if (d_changed)
   {
      // Apply changed Dirichlet BCs to T
      temp_u_gf.GetTrueDofs(u);

      // Here is where can set tmp_du_dt to be used if HO wanted
      // For higher-order, need T from n-1 timestep in addition to n timestep? is this worth doing?
      // Need to save previous timestep temperature in restarts

   }
   // Reassemble linear form b with updated coeffs + b_vec
   if (n_changed)
   {
      b->Assemble();
      b->ParallelAssemble(b_vec);
   }

}

void ConductionOperator::UpdateStiffness()
{    

   delete K_full;
   delete K;
   delete K_e;

   K_full = NULL;
   K = NULL;
   K_e = NULL;


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
   rhs.Add(1, b_vec);

   // Now have RHS without "eliminated" essential DOFs
   // ^this must be done depending on LHS
}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du/dt = M^{-1}(-Ku + boundary_terms)
   // for du_dt

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


ConductionOperator::~ConductionOperator()
{
   delete[] all_bdr_attr_markers;
   
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
