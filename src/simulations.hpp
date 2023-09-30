#pragma once

#include "mfem/mfem.hpp"
#include "option_structure.hpp"

class Simulation
{
    private:
    protected:
        mfem::Vector sol;
    public:
        virtual bool UsesMaterialProperty(MATERIAL_PROPERTY mat_prop) = 0;

        virtual void InitializeSolver() = 0;
        virtual bool Running() = 0;
        virtual void PreprocessIteration() = 0;
        virtual void Iterate() = 0;
        
        void SaveOldState();
        void ReloadOldState();
};

class UnsteadyHeatSimulation : public Simulation
{
    private:
    protected:
        double time;
        double dt;
        //ConductionOperator* oper;
    public:
        UnsteadyHeatSimulation() {};
        bool UsesMaterialProperty(MATERIAL_PROPERTY mat_prop);
        void InitializeSolver() {};
        bool Running() {};
        void PreprocessIteration() {};
        void Iterate() {};
        ~UnsteadyHeatSimulation() {};
};