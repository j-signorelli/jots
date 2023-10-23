#pragma once

#include "mfem/mfem.hpp"

#include "conduction_operator.hpp"
#include "option_structure.hpp"
#include "config_file.hpp"


class Simulation
{
    private:
    protected:
        const std::string solution_name;
        mfem::Vector u;
        mfem::Vector u_saved;

    public:
        Simulation(const std::string in_name, const int in_it_num, const mfem::ParGridFunction* u_0);
        std::string GetSolutionName() { return solution_name; };
        void SaveOldState() { u_saved = u; };
        void ReloadOldState() { u = u_saved;};
        void SetFromGF(const mfem::ParGridFunction* in_u) { in_u->GetTrueDofs(u); };
        const mfem::Vector& GetConstSolutionRef() const { return u; };
        
        virtual bool IsRunning() = 0;
        virtual void PreprocessIteration() = 0;
        virtual void Iterate() = 0;
        
};

class UnsteadyHeatSimulation : public Simulation
{
    private:
    protected:
        static const double TIME_TOLERANCE;

        const bool constant_stiffness;
        const bool constant_mass;
        const bool constant_neumann;
        const double tf;

        double& time;
        double& dt;

        mfem::ODESolver* ode_solver;
        ConductionOperator* oper;
    public:
        UnsteadyHeatSimulation(const mfem::ParGridFunction* u_0, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty* const* mat_props, ParFiniteElementSpace &f, double& in_time, double& dt);

        bool IsRunning();
        void PreprocessIteration();
        void Iterate();
        ~UnsteadyHeatSimulation();
};