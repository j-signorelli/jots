#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "boundary_condition.hpp"
#include "material_property.hpp"
#include "jots_iterator.hpp"
#include "helper_functions.hpp"
#include "jots_nlfis.hpp"
using namespace mfem;


// Reduced sytem operator for A(u,t) ---> See README.md
class ReducedSystemOperatorA : public Operator
{
    protected:
        const Operator* K; // not owned
        const ParNonlinearForm* B; // not owned - !NULL if rho(u) or C(u)
        const ParNonlinearForm* N; // not owned - !NULL if NL Neumann (rho(u) or C(u))
        const Vector* N_vec; // not owned - !NULL if linear Neummann
        mutable HypreParMatrix* Jacobian;
    public:
        ReducedSystemOperatorA(Operator* K_, ParNonlinearForm* B_, ParNonlinearForm* N_);  // For non-constant rhoC
        ReducedSystemOperatorA(Operator* K_, const Vector& N_vec_); // For constant rhoC
        void Mult(const Vector &u, Vector &y) const;
        Operator& GetGradient(const Vector &u) const;
        ~ReducedSystemOperatorA() { delete K; delete B; delete N; };
}

// Reduced system operator for R(u,t) ---> See README.md
class ReducedSystemOperatorR : public Operator
{
    private:
    protected:
        ParBilinearForm& m;
        ReducedSystemOperatorA& A;
        double& dt;
    public:
        ReducedSystemOperatorR(ParBilinearForm& m_, ReducedSystemOperatorA& A_, double& dt_) : Operator(), m(m_), A(A_), dt(dt_) {}; 
        void Mult(const Vector &k, Vector &y) const;
        Operator& GetGradient(const Vector &k) const;
}

class ConductionOperator : public TimeDependentOperator, public JOTSIterator
{
    protected:
        // Custom coefficients
        class AOverBCCoefficient : public Coefficient
        {
            private:
                Coefficient* A,B,C;
            public:
                AOverBCCoefficient(Coefficient& A_, Coefficient& B, Coefficient& C) : A(A_),  B(B_), C(C_) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);

        };
        class dAOverBCCoefficient : public Coefficient
        {
            private:
                Coefficient& A,dA,B,dB,C,dC;
            public:
                dAOverBCCoefficient(Coefficient& A_, Coefficient& dA_, Coefficient& B, Coefficient& dB_, Coefficient& C Coefficient& dC_) : A(A_), dA(dA_), B(B_), dB(dB_), C(C_), dC(dC_) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };

        double& time;
        double& dt;
        
        mfem::ODESolver* ode_solver;

        IterativeSolver *expl_solver;    // Solver for explicit time integration
        HypreSmoother expl_prec; // Preconditioner for the mass matrix M
        IterativeSolver *impl_solver;   // Solver for implicit time integration
        HypreSmoother impl_prec; // Preconditioner for the implicit solver
        JOTSNewtonSolver newton_solver;
        
        AOverBCCoefficient diffusivity, g_over_rhoC, one_over_rhoC;
        dAOverBCCoefficient d_diffusivity, dg_over_rhoC, done_over_rhoC;
        ProductCoefficient beta, d_beta;
        

        ParBilinearForm M;
        Operator* K;
        ParNonlinearForm* B;
        ParNonlinearForm* N;
        ReducedSystemOperatorA* A;

        HypreParMatrix *M_full; // full - no essential DOFs eliminated
        HypreParMatrix *K_full; // full - no essential DOFs eliminated

        HypreParMatrix *M; // ess DOFs eliminated - Operator for explicit time-integration
        HypreParMatrix *K; // ess DOFs eliminated
        HypreParMatrix *A; // ess DOFs eliminated - Operator for implicit time integration (A = M + dt K)
        
        HypreParMatrix *M_e; // eliminated part of M, M_full = M + M_e
        HypreParMatrix *K_e; // eliminated part of K, K_full = K + K_e
        HypreParMatrix *A_e; // eliminated part of A, A_full = A + A_e - required for setting BCs in implicit time integration
        mutable Vector rhs;
        mutable bool mass_updated;

        void CalculateRHS(const Vector &u) const;

        void ReassembleMass();

        void ReassembleStiffness();

    public:
        // Note: bdr attributes array cannot be constant. May move into BoundaryCondition class in future
        ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref);

        //--------------------------------------------------------------------------------------
        // TimeDependentOperator Function implementations:
        void Mult(const Vector &u, Vector &du_dt) const;

        /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
             This is the only requirement for high-order SDIRK implicit integration.*/
        void ImplicitSolve(const double dt, const Vector &u, Vector &k);

        //--------------------------------------------------------------------------------------
        // JOTSIterator Function implementations:
        void Iterate(mfem::Vector& u);

        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};

        ~ConductionOperator();
};