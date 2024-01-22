
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
        mfem::Array<MaterialProperty*> mps;
        const bool iterate_on_k;
        const mfem::Vector* u_n;
        const double* dt;
    public:
        JOTSNewtonSolver(MPI_Comm comm, const bool k) : mfem::NewtonSolver(comm), iterate_on_k(k), u_n(nullptr), dt(nullptr) {};
        void AddMaterialProperty(MaterialProperty& mp) { mps.Append(&mp); };
        void SetParameters(const mfem::Vector* u_n_, const double* dt_) { u_n = u_n_; dt = dt_; }
        void ProcessNewState(const mfem::Vector& x) override;
        void Mult(const mfem::Vector &b, mfem::Vector &x) const override;
}

#include "helper_functions.inl"