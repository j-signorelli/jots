#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "config_file.hpp"
#include "conduction_operator.hpp"
#include "conductivity_model.hpp"

class JOTSDriver
{
    protected:
        const std::string line = "-------------------------------------------------------------------------------------------";

	    int rank;
        int size;

        precice::SolverInterface* interface;

        Config* user_input;
        BoundaryCondition** boundary_conditions;
        ConductivityModel* cond_model;

        int dim;
        mfem::ODESolver* ode_solver;
        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* fespace;
        mfem::ParGridFunction* T_gf;
        mfem::Vector T;

        
        ConductionOperator* oper;

        double dt;

    public:
        JOTSDriver(const char* input_file, int myid, int num_procs);
        void Run();

        ~JOTSDriver();

    private:
};