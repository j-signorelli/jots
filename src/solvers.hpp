#pragma once

#include "mfem/mfem.hpp"


class Solver
{
    private:
    protected:
        mfem::Vector sol;
        bool initialized;
    public:
        Solver(): initialized(false) {}; // TODO: How to not be able to run any implemented functions if not initialized?
        virtual void InitializeSolver() = 0;
        virtual bool Running() = 0;
        virtual void PreprocessIteration() = 0;
        virtual void Iterate() = 0;
        
        void SaveOldState();
        void ReloadOldState();
}

class UnsteadyHeatSolver : public Solver
{
    private:
    protected:
        double time;
        double dt;
        ConductionOperator* oper;
    public:
        UnsteadyHeatSolver();
        void InitializeSolver();
        bool Running();
        void PreprocessIteration();
        void Iterate();
        ~UnsteadyHeatSolver();
}