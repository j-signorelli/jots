#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "boundary_condition.hpp"
#include "conductivity_model.hpp"

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
   const Config* user_input; // Not allocated here
   BoundaryCondition** boundary_conditions; // Not allocated here
   const bool constant_cond_model;

   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // list of essential true dofs
   Array<int>* all_bdr_attr_markers;

   IterativeSolver *expl_solver;    // Solver for explicit time integration
   HypreSmoother expl_prec; // Preconditioner for the mass matrix M
   IterativeSolver *impl_solver;   // Solver for implicit time integration
   HypreSmoother impl_prec; // Preconditioner for the implicit solver

   
   ParBilinearForm *m;
   ParBilinearForm *k;
   ParLinearForm *b;

   HypreParMatrix *M_full; // full - no essential DOFs eliminated
   HypreParMatrix *K_full; // full - no essential DOFs eliminated
   Vector b_vec; // full - Neumann BCs

   HypreParMatrix *M; // ess DOFs eliminated - Operator for explicit time-integration
   HypreParMatrix *K; // ess DOFs eliminated
   HypreParMatrix *A; // ess DOFs eliminated - Operator for implicit time integration (A = M + dt K)
   
   HypreParMatrix *M_e; // eliminated part of M, M_full = M + M_e
   HypreParMatrix *K_e; // eliminated part of K, K_full = K + K_e
   HypreParMatrix *A_e; // eliminated part of A, A_full = A + A_e - required for setting BCs in implicit time integration
   mutable Vector rhs; // = -KT + Neumann

   void PreprocessBCs();

   void PreprocessStiffness(const ConductivityModel* in_cond);
   
   void PreprocessSolver();

   // Apply the given boundary conditions
   void ApplyBCs(Vector &u);

   /// Update the diffusion BilinearForm K
   void UpdateStiffness();

   void CalculateRHS(const Vector &u) const;

public:
   ConductionOperator(const Config* in_config, BoundaryCondition** in_bcs, const ConductivityModel* in_cond, ParFiniteElementSpace &f, double t_0);

   void Mult(const Vector &u, Vector &du_dt) const;
   
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   void PreprocessIteration(Vector &u);
   
   ~ConductionOperator();
};