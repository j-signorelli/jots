#include "conduction_operator.hpp"

using namespace std;


/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
ConductionOperator::ConductionOperator(Config* in_config, ParFiniteElementSpace &f, const Vector &u, double t_0)
   :  TimeDependentOperator(f.GetTrueVSize(), t_0),
      fespace(f),
      M(NULL), 
      K(NULL),
      T(NULL),
      current_dt(t_0),
      M_solver(f.GetComm()),
      T_solver(f.GetComm()),
      user_input(in_config),
      z(height)
{
   const double rel_tol = 1e-8;


   // Set up parallel bilinear form
   M = new ParBilinearForm(&fespace);


   // Add the domain integral with included rho*Cp coefficient everywhere
   ConstantCoefficient rhoCp(user_input->GetDensity()*user_input->GetCp());
   M->AddDomainIntegrator(new MassIntegrator(rhoCp));
   M->Assemble(0); // keep sparsity pattern of M and K the same <-- TODO: understand why 0 is here
   
   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(user_input->GetBCCount());
   dbc_bdr = 0; // Initialize with all attributes set to non-essential = 0
   
   // Loop through vector boundary conditions, update array if needed (1=Dirichlet)
   for (size_t i = 0; i < user_input->GetBCCount(); i++)
      if (user_input->GetBCs()[i]->GetType() == BOUNDARY_CONDITION::ISOTHERMAL)
         dbc_bdr[i] = 1;
   
   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
   
   // Form the linear system matrix/operator accounting for the essential true DOFs, store in Mmat
   M->FormSystemMatrix(ess_tdof_list, Mmat);
   
   // Set up the solver we will use to invert Mmat
   M_solver.iterative_mode = false; // If true, would use second argument of Mult() as initial guess; here it is set to false
   M_solver.SetRelTol(rel_tol); // Sets "relative tolerance" of iterative solver
   M_solver.SetAbsTol(0.0); // Sets "absolute tolerance" of iterative solver
   M_solver.SetMaxIter(100); // Sets maximum number of iterations
   M_solver.SetPrintLevel(0); // Print all information about detected issues
   M_prec.SetType(HypreSmoother::Jacobi); // Set type of preconditioning (relaxation type) 
   M_solver.SetPreconditioner(M_prec); // Set preconditioner to matrix inversion solver
   M_solver.SetOperator(Mmat); // Set solver for given operator

   // Now set up the solver we will use for ____ (I have no clue yet)
   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec); // Use default preconditioning as none was set

}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep

   // Apply thermal conductivity model - cannot do! const
   //SetThermalConductivities(u);

   // Apply boundary conditions - cannot do! const
   //ApplyBCs(u);


   // ^Must do externally

   // Complete multiplication
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   M_solver.Mult(z, du_dt);
}


void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep

   // Apply thermal conductivity model - cannot do! const
   //SetThermalConductivities(u);

   // Apply boundary conditions - cannot do! const
   //ApplyBCs(u);


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

void ConductionOperator::SetThermalConductivities(const Vector &u)
{
   delete K;
   K = new ParBilinearForm(&fespace);

   // Apply thermal conductivity model, setup K appropriately
   Coefficient* u_coeff = user_input->GetConductivityModel()->ApplyModel(&fespace,u);

   K->AddDomainIntegrator(new DiffusionIntegrator(*u_coeff));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);
   delete T;
   T = NULL; // re-compute T on the next ImplicitSolve
}

void ApplyBCs(const Vector &u)
{
   //TODO:
}

ConductionOperator::~ConductionOperator()
{
   delete T;
   delete M;
   delete K;
}
