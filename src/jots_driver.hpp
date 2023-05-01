#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "config_file.hpp"
#include "conduction_operator.hpp"
#include "conductivity_model.hpp"
#include "precice_adapter.hpp"

class JOTSDriver
{
    protected:
        const std::string line = "-------------------------------------------------------------------------------";

	    const int rank;
        const int size;

        int dim;
        int it_num;
        double time;
        double dt;
        double tf;

        PreciceAdapter* adapter;

        Config* user_input;
        BoundaryCondition** boundary_conditions;
        ConductivityModel* cond_model;

        mfem::ODESolver* ode_solver;
        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* fespace;

        ConductionOperator* oper;

        mfem::ParGridFunction* T_gf;
        mfem::Vector T;


    public:
        JOTSDriver(const char* input_file, const int myid, const int num_procs);
        void Run();

        ~JOTSDriver();

    private:
};