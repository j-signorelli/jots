
#pragma once
#include "mfem/mfem.hpp"

#include "material_property.hpp"
#include "option_structure.hpp"

namespace Factory
{

mfem::IterativeSolver* GetSolver(std::string solver_label, MPI_Comm comm_);

mfem::HypreSmoother::Type GetPrec(std::string prec_label);

mfem::ODESolver* GetODESolver(std::string time_scheme_label);

}

template<typename Key, typename Value>
std::vector<Key> GetKeyVector(std::map<Key, Value> in_map);

class JOTSNewtonSolver : public mfem::NewtonSolver
{
    protected:
        mfem::Array<MaterialProperty*> mps; // None owned
        bool iterate_on_k;
    public:
        JOTSNewtonSolver(MPI_Comm comm) : mfem::NewtonSolver(comm), iterate_on_k(false) {};
        void AddMaterialProperty(MaterialProperty& mp) { mps.Append(&mp); };
        void SetOperator(const Operator& op) override;
        void ProcessNewState(const mfem::Vector& x) override;
}

// Use for Operators R(k) where k=dudt
class JOTS_k_Operator : public mfem::Operator
{
    protected:
        const mfem::Vector* u_n; // not owned
        const double* dt; // not owned
    public:
        JOTS_k_Operator(int s=0) : Operator(s), u_n(nullptr), dt(nullptr) {};
        JOTS_k_Operator(int h, int w) : Operator(h,w), u_n(nullptr), dt(nullptr) {}; 
        void SetParameters(const mfem::Vector* u_n_, const double* dt_) { u_n = u_n_; dt = dt_; }
        const Vector& Get_u_n() const { return *u_n; };
        const double& Get_dt() const {return *dt; };
}

#include "helper_functions.inl"