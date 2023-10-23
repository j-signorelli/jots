#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "config_file.hpp"
#include "material_property.hpp"
#include "precice_adapter.hpp"
#include "output_manager.hpp"
#include "simulations.hpp"

class JOTSDriver
{   
    private:

        void ProcessFiniteElementSetup();
        void ProcessMaterialProperties();
        void ProcessTimeIntegration();
        void ProcessPrecice();
        void ProcessBoundaryConditions();
        void PrintLinearSolverSettings();
        void PrintOutput();


    protected:
        static const std::string LINE;

	    const int rank;
        const int size;
        MPI_Comm comm;

        int dim;
        int it_num;
        double time;
        double precice_dt;
        double precice_saved_time;
        double precice_saved_it_num;
        double dt;
        double tf;

        Simulation* sim;

        PreciceAdapter* adapter;

        const Config& user_input;

        BoundaryCondition** boundary_conditions;
        Array<int>* all_bdr_attr_markers;
        bool initialized_bcs;

        MaterialProperty** mat_props;

        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* fespace;

        OutputManager* output;
        
        mutable mfem::ParGridFunction* temp_u_gf;

    public:
        JOTSDriver(const Config& input, const int myid, const int num_procs, MPI_Comm in_comm=MPI_COMM_WORLD);
        
        void UpdateMatProps();

        void UpdateBCs();

        void PreprocessIteration();

        void Iteration();

        void PostprocessIteration();

        void Run();
        
        const OutputManager* GetOutputManager() { return output; };

        ~JOTSDriver();
};