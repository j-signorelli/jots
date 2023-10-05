#pragma once

#include "mfem/mfem.hpp"
#include "option_structure.hpp"

class Simulation
{
    private:
    protected:
        const std::string solution_name;
        mfem::Vector u;
        mfem::Vector u_saved;
    public:
        std::string GetSolutionName() { return solution_name; };
        //mfem::Vector& GetSolutionVectorRef()
        void SaveOldState() { u_saved = u; };
        void ReloadOldState() { u = u_saved;};
        void SetFromGF(mfem::ParGridFunction* in_u) { in_u->GetTrueDofs(u); };

        virtual void InitializeSolver() = 0;
        virtual bool Running() = 0;
        virtual void PreprocessIteration() = 0;
        virtual void Iterate() = 0;
        
};

class UnsteadyHeatSimulation : public Simulation
{
    private:
    protected:
        double time;
        double dt;
        //ConductionOperator* oper;
    public:
        UnsteadyHeatSimulation() : solution_name("Temperature") {};

        void InitializeSolver() {};
        bool Running() {};
        void PreprocessIteration() {};
        void Iterate() {};
        ~UnsteadyHeatSimulation() {};
};