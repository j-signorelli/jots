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
ConductionOperator::ConductionOperator(Config* in_config, ParFiniteElementSpace &f, double t_0)
   :  TimeDependentOperator(f.GetTrueVSize(), t_0),
      fespace(f),
      k_coeff(NULL),
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
      user_input(in_config)
{

   PreprocessBCs();
   
   PreprocessStiffness();

   PreprocessSolver();

}

void ConductionOperator::PreprocessBCs()
{

   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(user_input->GetBCCount());
   dbc_bdr = 0; // Initialize with all attributes set to non-essential = 0

   // Further, create bdr_attr_marker arrays and FunctionCoefficient arrays for all BCs
   all_bdr_attr_markers = new Array<int>[user_input->GetBCCount()];
   all_bdr_coeffs = new Coefficient*[user_input->GetBCCount()];
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      Array<int> bdr_attr(user_input->GetBCCount());
      bdr_attr = 0;
      bdr_attr[i] = 1;
      all_bdr_attr_markers[i] = bdr_attr;

      all_bdr_coeffs[i] = NULL;

      // Update DBC list if essential
      if (user_input->GetBCs()[i]->IsEssential())
         dbc_bdr[i] = 1;
   }

   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // Create linear form for Neumann BCs
   b = new ParLinearForm(&fespace);

   // Loop through all BCs
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {  
      // Get the coefficients for all BCs
      all_bdr_coeffs[i] = user_input->GetBCs()[i]->GetCoefficient();
   
      // Add Neumann BCs to linear form b
      if (! user_input->GetBCs()[i]->IsEssential())
      {
         b->AddBoundaryIntegrator(new BoundaryLFIntegrator(*all_bdr_coeffs[i]), all_bdr_attr_markers[i]);// Add boundary integrator to boundary
      }
   }


   // Assemble b - TODO: Verify that time-dependent coeffs don't crash if SetTime not called first
   b->Assemble();

   // Set b_vec
   b->ParallelAssemble(b_vec);
}

void ConductionOperator::PreprocessStiffness()
{  

   // Assemble the parallel bilinear form for stiffness matrix
   k = new ParBilinearForm(&fespace);

   // Get thermal conductivity coefficient
   k_coeff = user_input->GetConductivityModel()->GetCoefficient();
   
   // Add domain integrator to the bilinear form
   k->AddDomainIntegrator(new DiffusionIntegrator(*k_coeff));


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
   
   // Create stiffness matrix w/ removed essential DOFs + remove from M
   M_e = k->ParallelEliminateTDofs(ess_tdof_list, *M);

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

void ConductionOperator::PreprocessIteration(Vector &u, double curr_time)
{
   // Apply BCs
   ApplyBCs(u, curr_time);

   //Calculate thermal conductivities
   SetThermalConductivities(u, curr_time);

}

void ConductionOperator::ApplyBCs(Vector &u, double curr_time)
{

   // Update time-dependent coefficients with current time
   // Also project essential coefficients onto u
   ParGridFunction temp_u_gf(&fespace);
   temp_u_gf.SetFromTrueDofs(u);

   bool n_changed = false;
   bool d_changed = false;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      if (!user_input->GetBCs()[i]->IsConstant()) // If not constant in time
      {
         // Set coefficient's time to current time
         all_bdr_coeffs[i]->SetTime(curr_time);
         n_changed = true;
      }

      if (user_input->GetBCs()[i]->IsEssential())
      {
         // Project correct values on boundary for essential BCs
         d_changed = true;
         temp_u_gf.ProjectBdrCoefficient(*all_bdr_coeffs[i], all_bdr_attr_markers[i]);
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

void ConductionOperator::SetThermalConductivities(const Vector &u, double curr_time)
{    

   // Update matrix K IF TIME-DEPENDENT k Or if not yet instantiated
   if (!K || !user_input->GetConductivityModel()->IsConstant())
   {  


      delete K_full;
      delete K;
      delete K_e;


      k_coeff->SetTime(curr_time);     
      k->Update(); // delete old data (M and M_e)
      k->Assemble(0);
      k->Finalize(0);


      // Create full stiffness matrix w/o removed essential DOFs
      K_full = k->ParallelAssemble();
      K = new HypreParMatrix(*K_full);
   
      // Create stiffness matrix w/ removed essential DOFs
      K_e = k->ParallelEliminateTDofs(ess_tdof_list, *K);

   }
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

   // Calculate RHS pre-essential update
   CalculateRHS(u);

   // Calculate LHS w/o essential DOFs eliminated at first
   A = Add(1.0, *M, dt, *K);

   impl_solver->SetOperator(*A);

   // Calculate LHS eliminated part using already calculated M_e's from the individual bilinear forms
   A_e = Add(1.0, *M_e, dt, *K_e);

   du_dt.SetSubVector(ess_tdof_list, 0.0); // **Assume essential DOFs du_dt = 0.0

   // Eliminate BCs from RHS given LHS
   A->EliminateBC(*A_e, ess_tdof_list, du_dt, rhs);

   // Solve for du_dt
   impl_solver->Mult(rhs, du_dt);
}


ConductionOperator::~ConductionOperator()
{
   //delete T;
   delete[] all_bdr_attr_markers;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
      delete all_bdr_coeffs[i];
   delete[] all_bdr_coeffs;
   delete k_coeff;
   
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
