#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "config_file.hpp"
#include "conduction_operator.hpp"
#include "material_property.hpp"
#include "precice_adapter.hpp"
#include "output_manager.hpp"

class JOTSDriver
{
    protected:
        static const std::string LINE;
        static const double TIME_TOLERANCE;

	    const int rank;
        const int size;

        int dim;
        int it_num;
        double time;
        double dt;
        double tf;

        PreciceAdapter* adapter;

        Config& user_input;

        BoundaryCondition** boundary_conditions;
        Array<int>* all_bdr_attr_markers;
        bool initialized_bcs;
        
        MaterialProperty* k_prop;
        MaterialProperty* C_prop;

        mfem::ODESolver* ode_solver;
        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* fespace;

        ConductionOperator* oper;
        OutputManager* output;

        mfem::Vector T;
        
        mutable mfem::ParGridFunction* temp_T_gf;

        void UpdateMatProps();

        void PreprocessIteration();

    public:
        JOTSDriver(const char* input_file, const int myid, const int num_procs);
        void Run();


        ~JOTSDriver();

    private:
};