#pragma once
#include <cstdio>

#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "config_file.hpp"
#include "material_property.hpp"
#include "jots_precice.hpp"
#include "output_manager.hpp"
#include "jots_iterator.hpp"
#include "linear_conduction_operator.hpp"
#include "nl_conduction_operator.hpp"
#include "steady_conduction_operator.hpp"
#include "boundary_condition.hpp"

class JOTSDriver
{   
    private:

        void ProcessFiniteElementSetup();
        void ProcessMaterialProperties();
        void ProcessTimeIntegration();
        void ProcessPrecice();
        void ProcessBoundaryConditions();
        void PrintLinearSolverSettings();
        void PrintNewtonSolverSettings();
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
        int max_timesteps;

        JOTSIterator* jots_iterator[PHYSICS_TYPE_SIZE];
        mfem::Vector* u[PHYSICS_TYPE_SIZE];

        JOTSSolverInterface* precice_interface;

        const Config& user_input;

        BoundaryCondition** boundary_conditions[PHYSICS_TYPE_SIZE];
        Array<int>* all_bdr_attr_markers;
        bool initialized_bcs;

        MaterialProperty* mat_props[MATERIAL_PROPERTY_SIZE];

        mfem::ParMesh* pmesh;
        mfem::FiniteElementCollection* fe_coll;
        mfem::ParFiniteElementSpace* scalar_fespace;
        mfem::ParFiniteElementSpace* vector_fespace;

        // Array associating each physics-type to its FESpace
        // Pointers to either scalar_fespace or vector_fespace, or nullptr if not solving
        mfem::ParFiniteElementSpace* fespace[PHYSICS_TYPE_SIZE];

        OutputManager* output;

        
        mutable mfem::ParGridFunction* u_0_gf[PHYSICS_TYPE_SIZE];

    public:
        JOTSDriver(const Config& input, const int myid, const int num_procs, MPI_Comm in_comm=MPI_COMM_WORLD);
        
        void UpdateMatProps(const bool apply_changes);

        void UpdateAndApplyBCs();

        void PreprocessIteration();

        void Iteration();

        void PostprocessIteration();

        void Run();
        
        OutputManager* GetOutputManager() { return output; };

        ~JOTSDriver();
};