#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"

#include "config_file.hpp"
#include "conduction_operator.hpp"
#include "conductivity_model.hpp"
using namespace mfem;

class JOTSDriver
{
    protected:
        const std::string line = "-------------------------------------------------------------------------------------------";

	    int rank;

        Config* user_input;
        int dim;
        ODESolver* ode_solver;
        ParMesh* pmesh;
        FiniteElementCollection* fe_coll;
        ParFiniteElementSpace* fespace;
        ParGridFunction* T_gf;
        Vector T;

        
        ConductionOperator* oper;


    public:
        JOTSDriver(const char* input_file, int myid);
        void Run();

        ~JOTSDriver();

    private:
};