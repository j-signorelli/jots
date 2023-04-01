#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"

using namespace mfem;

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
class ConductionOperator : public TimeDependentOperator
{
protected:
   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // list of essential true dofs

   ParBilinearForm *M;
   ParBilinearForm *K;
   ParLinearForm *b; // Linear form for Neumann BCs


   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   //HypreParMatrix *T; // T = M + dt K
   double current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   //CGSolver T_solver;    // Implicit solver for T = M + dt K
   //HypreSmoother T_prec; // Preconditioner for the implicit solver

   Config* user_input; // Not allocated here

   mutable Vector z; // auxiliary vector
   
public:
   ConductionOperator(Config* in_config, ParFiniteElementSpace &f, double t_0);

   void Mult(const Vector &u, Vector &du_dt) const;
   
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   //void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   /// Update the diffusion BilinearForm K using the given true-dof vector `u` based on specified model.
   void SetThermalConductivities(const Vector &u);

   // Apply the given boundary conditions
   void ApplyBCs(ParGridFunction* u_gf, Vector &u);

   ~ConductionOperator();
};