#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "boundary_condition.hpp"
#include "material_property.hpp"
#include "jots_iterator.hpp"
#include "helper_functions.hpp"

using namespace mfem;


class ConductionOperator : public TimeDependentOperator, public JOTSIterator
{
protected:

   const double& tf;
   double& time;
   double& dt;

   mfem::ProductCoefficient rho_C;

   ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // list of essential true dofs

   mfem::ODESolver* ode_solver;

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
   mutable bool mass_updated;

   void PreprocessBCs(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers);

   void PreprocessMass();

   void PreprocessStiffness(const MaterialProperty* k_prop);
   
   void PreprocessSolver(const Config& in_config);

   void CalculateRHS(const Vector &u) const;

   void ReassembleMass();
   
   void ReassembleStiffness();

public:
   // Note: bdr attributes array cannot be constant. May move into BoundaryCondition class in future
   ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty* rho_prop, const MaterialProperty* C_prop, const MaterialProperty* k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref, const double& tf_ref);
   
   //--------------------------------------------------------------------------------------
   // TimeDependentOperator Function implementations:
   void Mult(const Vector &u, Vector &du_dt) const;
   
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   //--------------------------------------------------------------------------------------
   // JOTSIterator Function implementations:
   bool IsNotComplete() const;

   void Iterate(mfem::Vector& u);

   void ProcessMatPropUpdate(MATERIAL_PROPERTY mp);

   void UpdateNeumann();

   ~ConductionOperator();
};