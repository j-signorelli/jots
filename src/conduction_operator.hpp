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
};

// Reduced system operator for R(k) ---> See README.md
class ReducedSystemOperatorR : public JOTS_k_Operator // Use JOTS_k_Operator
{
    private:
    protected:
        const HypreParMatrix& M_mat;
        const ReducedSystemOperatorA& A;
        const Array<int>& ess_tdof_list;
        mutable HypreParMatrix* Jacobian; // owned
        mutable Vector z;
    public:
        ReducedSystemOperatorR(const HypreParMatrix& M_mat_, const ReducedSystemOperatorA& A_, const Array<int>& ess_tdof_list_) : JOTS_k_Operator(M_mat_.Height()), M_mat(M_mat_), A(A_), ess_tdof_list(ess_tdof_list_), Jacobian(nullptr), z(height) {}; 
        void Mult(const Vector &k, Vector &y) const;
        Operator& GetGradient(const Vector &k) const;
        ~ReducedSystemOperatorR() { delete Jacobian; };
};

class ConductionOperator : public TimeDependentOperator, public JOTSIterator
{
    protected:
        class DiffusivityCoefficient : public Coefficient
        {
            protected:
                Coefficient &k, &rho, &C;
            public:
                DiffusivityCoefficient(const MaterialProperty& k_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : k(k_.GetCoeffRef()), rho(rho_.GetCoeffRef()), C(C_.GetCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        class dDiffusivityCoefficient : public Coefficient
        {
            protected:
                Coefficient &k, &dk, &rho, &drho, &C, &dC;
            public:
                dDiffusivityCoefficient(const MaterialProperty& k_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : k(k_.GetCoeffRef()), dk(k_.GetDCoeffRef()), rho(rho_.GetCoeffRef()), drho(rho_.GetDCoeffRef()), C(C_.GetCoeffRef()), dC(C_.GetDCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        class NeumannCoefficient : public Coefficient
        {
            protected:
                Coefficient &g, &rho, &C;
            public:
                NeumannCoefficient(Coefficient& g_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : g(g_), rho(rho_.GetCoeffRef()), C(C_.GetCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        class dNeumannCoefficient : public Coefficient
        {   
            protected:
                Coefficient &g, &rho, &drho, &C, &dC;
            public:
                dNeumannCoefficient(Coefficient& g_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : g(g_), rho(rho_.GetCoeffRef()), drho(rho_.GetDCoeffRef()), C(C_.GetCoeffRef()), dC(C_.GetDCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        class BetaCoefficient : public Coefficient
        {   
            protected:
                Coefficient &k, &dk, &rho, &drho, &C, &dC;
            public:
                BetaCoefficient(const MaterialProperty& k_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : k(k_.GetCoeffRef()), dk(k_.GetDCoeffRef()), rho(rho_.GetCoeffRef()), drho(rho_.GetDCoeffRef()), C(C_.GetCoeffRef()), dC(C_.GetDCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        class dBetaCoefficient : public Coefficient
        {
            protected:
                Coefficient &k, &dk, &rho, &drho, &d2rho, &C, &dC, &d2C;
            public:
                dBetaCoefficient(const MaterialProperty& k_, const MaterialProperty& rho_, const MaterialProperty& C_)
                : k(k_.GetCoeffRef()), dk(k_.GetDCoeffRef()), rho(rho_.GetCoeffRef()), drho(rho_.GetDCoeffRef()), d2rho(rho_.GetD2CoeffRef()), C(C_.GetCoeffRef()), dC(C_.GetDCoeffRef()), d2C(C_.GetD2CoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };

        const double &dt;
        mfem::ODESolver* ode_solver;

        IterativeSolver *lin_solver;    // Linear solver
        HypreSmoother lin_prec; // Preconditioner for linear solver
        JOTSNewtonSolver newton;
        
        DiffusivityCoefficient diffusivity;
        dDiffusivityCoefficient d_diffusivity;

        NeumannCoefficient g_over_rhoC;
        dNeumannCoefficient dg_over_rhoC;

        BetaCoefficient beta;
        dBetaCoefficient d_beta;

        ParBilinearForm M;
        Operator* K;
        ParNonlinearForm B;
        ParNonlinearForm N;
        ReducedSystemOperatorA* A;
        JOTS_k_Operator* R;
        HypreParMatrix M_mat;

        mutable Vector rhs;

        public:
            // Note: bdr attributes array cannot be constant. May move into BoundaryCondition class in future
            ConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& rho_prop, MaterialProperty& C_prop, MaterialProperty& k_prop, ParFiniteElementSpace &f, double &time_, double &dt_);

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