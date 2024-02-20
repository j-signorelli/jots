#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "thermal_boundary_condition.hpp"
#include "material_property.hpp"
#include "jots_iterator.hpp"
#include "jots_common.hpp"

using namespace mfem;


class LinearConductionOperator : public TimeDependentOperator, public JOTSIterator
{
protected:

	double& time;
	double& dt;

	mfem::ProductCoefficient rho_C;

	mfem::ODESolver* ode_solver;

	IterativeSolver *expl_solver;    // Solver for explicit time integration
	HypreSmoother expl_prec; // Preconditioner for the mass matrix M
	IterativeSolver *impl_solver;   // Solver for implicit time integration
	HypreSmoother impl_prec; // Preconditioner for the implicit solver

	
	ParBilinearForm M;
	ParBilinearForm K;

	HypreParMatrix *M_mat;
	HypreParMatrix *K_mat;
	HypreParMatrix *A_mat; // Operator for implicit time integration (A = M + dt K)
	
	mutable Vector rhs; // used for = -Ku + N

	void CalculateRHS(const Vector &u, Vector &y) const;

	void ReassembleM();

	void ReassembleK();
    
    void ReassembleA();

public:
	// Note: bdr attributes array cannot be constant. May move into BoundaryCondition class in future
	LinearConductionOperator(const Config &in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty &rho_prop, const MaterialProperty &C_prop, const MaterialProperty &k_prop, ParFiniteElementSpace &f, double &t_ref, double &dt_ref);

	//--------------------------------------------------------------------------------------
	// TimeDependentOperator Function implementations:
	void Mult(const Vector &u, Vector &du_dt) const;

	/** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
		 This is the only requirement for high-order SDIRK implicit integration.*/
	void ImplicitSolve(const double dt, const Vector &u, Vector &k);

	//--------------------------------------------------------------------------------------
	// JOTSIterator Function implementations:

	void Iterate(mfem::Vector& u);

	void ProcessMatPropUpdate(MATERIAL_PROPERTY mp);


	~LinearConductionOperator();
};