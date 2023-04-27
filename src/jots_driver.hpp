#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "config_file.hpp"
#include "conduction_operator.hpp"
#include "conductivity_model.hpp"
#include "solver_state.hpp"

class JOTSDriver
{
    protected:
        const std::string line = "-------------------------------------------------------------------------------------------";

	    const int rank;
        const int size;

        int dim;

        precice::SolverInterface* interface;

        Config* user_input;
        BoundaryCondition** boundary_conditions;
        ConductivityModel* cond_model;
        SolverState* state;

        mfem::ODESolver* ode_solver;
        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* fespace;

        ConductionOperator* oper;

    public:
        JOTSDriver(const char* input_file, const int myid, const int num_procs);
        void Run();

        ~JOTSDriver();

    private:
};