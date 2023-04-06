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
   Array<int>* all_bdr_attr_markers;
   Coefficient** all_bdr_coeffs;
   Coefficient* k_coeff;
   
   ParBilinearForm *m;
   ParBilinearForm *k;
   ParLinearForm *b;

   HypreParMatrix M;
   HypreParMatrix K;
   //HypreParMatrix *T; // T = M + dt K
   //double current_dt;

   FGMRESSolver solver;    // FMGRES solver for inverting the mass matrix M
   HypreSmoother prec; // Preconditioner for the mass matrix M

   //CGSolver T_solver;    // Implicit solver for T = M + dt K
   //HypreSmoother T_prec; // Preconditioner for the implicit solver

   Config* user_input; // Not allocated here

   //Vector* b_vec; // Vector for enforcing Neumann BCs
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
   void ApplyBCs(Vector &u);

   ~ConductionOperator();
};