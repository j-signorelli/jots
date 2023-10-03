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
#include "simulations.hpp"

class JOTSDriver
{   
    private:

        void ProcessFiniteElementSetup();
        void ProcessMaterialProperties();
        void ProcessTimeIntegration();
        void ProcessPrecice();
        void ProcessBoundaryConditions();
        void ProcessLinearSolverSettings();
        void ProcessOutput();

        //void InitializeSolver();


    protected:
        static const std::string LINE;
        static const double TIME_TOLERANCE;

        const SIMULATION_TYPE sim_type;
	    const int rank;
        const int size;
        MPI_Comm comm;

        int dim;
        int it_num;
        double time;
        double dt;
        double tf;

        Simulation* sim;

        PreciceAdapter* adapter;

        const Config& user_input;

        BoundaryCondition** boundary_conditions;
        Array<int>* all_bdr_attr_markers;
        bool initialized_bcs;

        std::map<MATERIAL_PROPERTY, MaterialProperty*> mat_props;

        mfem::ODESolver* ode_solver;
        mfem::IterativeSolver* lin_solver;
        mfem::HypreSmoother::Type prec;
        
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
        JOTSDriver(const Config& input, const int myid, const int num_procs, MPI_Comm in_comm=MPI_COMM_WORLD);
        void Run();
        
        OutputManager* GetOutputManager() { return output; };

        ~JOTSDriver();
};