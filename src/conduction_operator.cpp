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
      user_input(in_config),
      z(height)
{
   const double rel_tol = 1e-16;
   const double abs_tol = 1e-10;

   // Assemble parallel bilinear form for mass matrix now (this does not change in time)
   m = new ParBilinearForm(&fespace);

   // Add the domain integral with included rho*Cp coefficient everywhere
   ConstantCoefficient rhoCp(user_input->GetDensity()*user_input->GetCp());
   m->AddDomainIntegrator(new MassIntegrator(rhoCp));
   m->Assemble(0); // keep sparsity pattern of m and k the same
   
   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(in_config->GetBCCount());
   dbc_bdr = 0; // Initialize with all attributes set to non-essential = 0
   
   // Loop through vector boundary conditions, update array if needed (1=Dirichlet)
   // Further, create bdr_attr_marker arrays and Coefficient arrays for all BCs
   // ^This is required as everything in MFEM refers to references, so we must keep everything intact until simulation is complete
   all_bdr_attr_markers = new Array<int>[user_input->GetBCCount()];
   all_bdr_coeffs = new Coefficient*[user_input->GetBCCount()];
   
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {
      Array<int> bdr_attr(user_input->GetBCCount());
      bdr_attr = 0;
      bdr_attr[i] = 1;
      all_bdr_attr_markers[i] = bdr_attr;

      all_bdr_coeffs[i] = NULL;

      if (user_input->GetBCs()[i]->IsEssential())
         dbc_bdr[i] = 1;
   }


   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
   

   // Form the linear system matrix/operator, store in Mmat
   m->FormSystemMatrix(ess_tdof_list, M);
   
   // Set up the solver
   solver.iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false

   solver.SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
   solver.SetAbsTol(abs_tol); // Sets "absolute tolerance" of iterative solver
   solver.SetMaxIter(100); // Sets maximum number of iterations

   solver.SetPrintLevel(0); // Print all information about detected issues

   prec.SetType(HypreSmoother::Chebyshev); // Set type of preconditioning (relaxation type) 
   solver.SetPreconditioner(prec); // Set preconditioner to matrix inversion solver

   solver.SetOperator(M); // Set matrix/operator M as the LHS for the solver


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

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du/dt = M^{-1}(-Ku + boundary_terms)
   // for du_dt

   // Complete multiplication of Ku
   K.Mult(u, z);
   z.Neg();

   // Add boundary term (for Neumann BCs)
   //z += *b;
   //m->FormLinearSystem(ess_tdof_list, z, )
   // Solver M^-1, then multiply M^-1 * z
   solver.Mult(z, du_dt);

   //test:
   //du_dt.SetSubVector(ess_tdof_list, 0.0);
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

void ConductionOperator::SetThermalConductivities(const Vector &u)
{  

   delete K;
   delete k_coeff;

   // Assemble the parallel bilinear form for stiffness matrix
   k = new ParBilinearForm(&fespace);

   // Apply thermal conductivity model to get appropriate coefficient \lambda in Diffusion integrator
   k_coeff = user_input->GetConductivityModel()->ApplyModel(&fespace,u);

   // Add domain integrator to the bilinear form, and then obtain a system matrix K
   k->AddDomainIntegrator(new DiffusionIntegrator(*k_coeff));
   k->Assemble(0); // keep sparsity pattern of m and k the same
   k->FormSystemMatrix(ess_tdof_list, K); // Create Operator K

   //delete T;
   //T = NULL; // re-compute T on the next ImplicitSolve
}

void ConductionOperator::ApplyBCs(Vector &u)
{
   // Delete previous stuff
   delete b;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
      delete all_bdr_coeffs[i];
   
   // Create temp ParGridFunction u_gf with values from (updated) u
   ParGridFunction u_gf_temp(&fespace);
   u_gf_temp.SetFromTrueDofs(u);

   // Create linear form for Neumann BCs
   b = new ParLinearForm(&fespace);

   // Create ParGridFunction using 
   // Loop through all BCs
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
   {  
      // Get the coefficient of interest
      all_bdr_coeffs[i] = user_input->GetBCs()[i]->GetCoefficient();

      // Set appropriately
      if (user_input->GetBCs()[i]->IsEssential())
      {
         u_gf_temp.ProjectBdrCoefficient(*all_bdr_coeffs[i], all_bdr_attr_markers[i]);  // Update u
      } else
      {
         b->AddBoundaryIntegrator(new BoundaryLFIntegrator(*all_bdr_coeffs[i]), all_bdr_attr_markers[i]);// Add boundary integrator to boundary
      }


   }
   
   // Update vector u Dirichlet BC values
   u_gf_temp.GetTrueDofs(u);

   // Assemble the new linear form for Neumann BCs
   b->Assemble();


}

ConductionOperator::~ConductionOperator()
{
   //delete T;
   delete[] all_bdr_attr_markers;
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
      delete all_bdr_coeffs[i];
   delete k_coeff;
   delete[] all_bdr_coeffs;
   delete m;
   delete K;
   delete b;
}
