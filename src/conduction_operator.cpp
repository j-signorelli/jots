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
      m(NULL), 
      k(NULL),
      b(NULL),
      //T(NULL),
      //current_dt(t_0),
      solver(f.GetComm()),
      //T_solver(f.GetComm()),
      user_input(in_config)
{

   
   PreprocessBCs();
   
   PreprocessStiffness();

   PreprocessSolver();

   // Now set up the solver we will use for implicit time integration
   // TODO: not yet implemented (understood)
   /*
   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec); // Use default preconditioning as none was set
   */
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
}

void ConductionOperator::PreprocessStiffness()
{
   // Assemble the parallel bilinear form for stiffness matrix
   k = new ParBilinearForm(&fespace);

   // Get thermal conductivity coefficient
   k_coeff = user_input->GetConductivityModel()->GetCoefficient();
   
   // Add domain integrator to the bilinear form
   k->AddDomainIntegrator(new DiffusionIntegrator(*k_coeff));
   k->Assemble(0); // keep sparsity pattern of m and k the same
   k->Finalize(0);

}

void ConductionOperator::PreprocessSolver()
{  
   const double rel_tol = 1e-16;
   const double abs_tol = 1e-10;


   delete m;

   // Assemble parallel bilinear form for mass matrix
   m = new ParBilinearForm(&fespace);

   // Add the domain integral with included rho*Cp coefficient everywhere
   ConstantCoefficient rhoCp(user_input->GetDensity()*user_input->GetCp());
   m->AddDomainIntegrator(new MassIntegrator(rhoCp));
   m->Assemble(0); // keep sparsity pattern of m and k the same
   

   // Form the linear system matrix/operator, store in M
   m->FormSystemMatrix(ess_tdof_list, M);
   
   // Set up the solver
   solver.iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false

   solver.SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
   solver.SetAbsTol(abs_tol); // Sets "absolute tolerance" of iterative solver
   solver.SetMaxIter(100); // Sets maximum number of iterations

   solver.SetPrintLevel(0); // Print all information about detected issues

   prec.SetType(HypreSmoother::Chebyshev); // Set type of preconditioning (relaxation type) 
   solver.SetPreconditioner(prec); // Set preconditioner to matrix inversion solver

   solver.SetOperator(M); // Set matrix/operator M as the LHS for the solver - this current configuration does not allow this to vary in time

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

   bool changed = false;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      if (!user_input->GetBCs()[i]->IsConstant()) // If not constant in time
      {
         // Set coefficient's time to current time
         all_bdr_coeffs[i]->SetTime(curr_time);
         changed = true;
      }

      if (user_input->GetBCs()[i]->IsEssential())
      {
         // Project correct values on boundary for essential BCs
         //cout << "Projecting BDR coefficient " << all_bdr_coeffs[i]->GetValue()
         temp_u_gf.ProjectBdrCoefficient(*all_bdr_coeffs[i], all_bdr_attr_markers[i]);
      }
   }

   // Apply Dirichlet BCs to T
   temp_u_gf.GetTrueDofs(u);

   // Reassemble linear form b with updated coeffs
   if (changed)
      b->Assemble();

}

void ConductionOperator::SetThermalConductivities(const Vector &u, double curr_time)
{  
   // Update matrix K IF TIME-DEPENDENT k
   if (!user_input->GetConductivityModel()->IsConstant())
   {
      // TODO: May actually need to create new coefficient
      // Can I just change ptr of k_coeff?
      
      k_coeff->SetTime(curr_time);     
      k->Update();
      k->Assemble(0);
      k->Finalize(0);
   }
   //delete T;
   //T = NULL; // re-compute T on the next ImplicitSolve
}



void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du/dt = M^{-1}(-Ku + boundary_terms)
   // for du_dt

   // Create auxiliary ParLinearForm to store dual vector values
   ParLinearForm z(&fespace);

   // Complete multiplication of ParBilinearForm with primal vector for dual vector z
   k->Mult(u, z);

   z.Neg(); // z = -z

   // Add Neumann ParLinearForm term
   z.Add(1, *b);

   // Now have RHS (as required dual vector / ParLinearForm)!

   // Enforce appropriate du_dt=0 for Dirichlet BCs
   ParGridFunction tmp_du_dt(&fespace);
   tmp_du_dt = 0.0; // TODO: Actually calculate w/ backward differencing for unsteady Dirichlet

   // Auxiliary variables
   OperatorPtr A;
   Vector B, X;

   // Form Linear System
   m->FormLinearSystem(ess_tdof_list, tmp_du_dt, z, A, X, B);


   // Solver M^-1, then multiply M^-1 * z
   solver.Mult(B, X);

   
   //----------------------------------------------------------
   // TODO: Check this more below, Rob just does above Mult with X=du_dt
   //       I believe that below is more generalized
   // Recover the solution as a ParGridFunction
   m->RecoverFEMSolution(X, z, tmp_du_dt);

   // Set new DOFs
   tmp_du_dt.GetTrueDofs(du_dt);
   

}

/* TODO:
void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep

   // ^Must do externally
   if (!T)
   {
      T = Add(1.0, Mmat, dt, Kmat);
      current_dt = dt;
      T_solver.SetOperator(*T);
   }
   MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   T_solver.Mult(z, du_dt);
}
*/

ConductionOperator::~ConductionOperator()
{
   //delete T;
   delete[] all_bdr_attr_markers;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
      delete all_bdr_coeffs[i];
   delete k_coeff;
   delete[] all_bdr_coeffs;
   delete m;
   //delete K;
   delete b;
}
