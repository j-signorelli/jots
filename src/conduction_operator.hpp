#pragma once
#include <vector>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "boundary_condition.hpp"
#include "material_property.hpp"
#include "jots_iterator.hpp"
#include "jots_common.hpp"
#include "jots_nlfis.hpp"
using namespace mfem;


// Reduced sytem operator for A(u,t) ---> See README.md
class ReducedSystemOperatorA : public Operator
{
    protected:
        const Operator* K; // not owned - either a HypreParMatrix or ParNonlinearForm
        const ParNonlinearForm* B; // not owned - !NULL if rho(u) or C(u)
        const ParNonlinearForm* N; // not owned - !NULL if NL Neumann (rho(u) or C(u))
        const Vector* N_vec; // not owned - !NULL if linear Neummann
        mutable HypreParMatrix* Jacobian; // owned
    public:
        ReducedSystemOperatorA(const Operator* K_, const ParNonlinearForm* B_, const ParNonlinearForm* N_);  // For non-constant rhoC
        ReducedSystemOperatorA(const Operator* K_, const Vector* N_vec_); // For constant rhoC
        void Mult(const Vector &u, Vector &y) const;
        Operator& GetGradient(const Vector &u) const;
        ~ReducedSystemOperatorA() { delete Jacobian; };
}

// Reduced system operator for R(k) ---> See README.md
class ReducedSystemOperatorR : public JOTS_k_Operator // Use JOTS_k_Operator
{
    private:
    protected:
        const ParBilinearForm& M;
        const ReducedSystemOperatorA& A;
        const Array<int>& ess_tdof_list;
        mutable HypreParMatrix* Jacobian; // owned
        mutable Vector z;
    public:
        ReducedSystemOperatorR(ParBilinearForm& m_, ReducedSystemOperatorA& A_, const Array<int>& ess_tdof_list_) : JOTS_k_Operator(m_ParFESpace()->TrueVSize()), m(m_), A(A_), ess_tdof_list_(ess_tdof_list_), Jacobian(nullptr), z(height) {}; 
        void Mult(const Vector &k, Vector &y) const;
        Operator& GetGradient(const Vector &k) const;
        ~ReducedSystemOperatorR() { delete Jacobian; };
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
        
        mfem::ODESolver* ode_solver;

        IterativeSolver *lin_solver;    // Linear solver
        HypreSmoother lin_prec; // Preconditioner for linear solver
        JOTSNewtonSolver newton;
        
        AOverBCCoefficient diffusivity, g_over_rhoC, one_over_rhoC;
        dAOverBCCoefficient d_diffusivity, dg_over_rhoC, done_over_rhoC;
        ProductCoefficient beta, d_beta;

        ParBilinearForm M;
        Operator* K;
        ParNonlinearForm* B;
        ParNonlinearForm* N;
        ReducedSystemOperatorA* A;
        JOTS_k_Operator* R;
        HypreParMatrix M_mat;

        mutable Vector rhs;

        public:
            // Note: bdr attributes array cannot be constant. May move into BoundaryCondition class in future
            ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double& t_ref, double& dt_ref);

            //--------------------------------------------------------------------------------------
            // TimeDependentOperator Function implementations:
            void Mult(const Vector &u, Vector &k) const;

            /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
                 This is the only requirement for high-order SDIRK implicit integration.*/
            void ImplicitSolve(const double dt, const Vector &u, Vector &k);

            //--------------------------------------------------------------------------------------
            // JOTSIterator Function implementations:
            void Iterate(mfem::Vector& u);

            void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};

        ~ConductionOperator();
};