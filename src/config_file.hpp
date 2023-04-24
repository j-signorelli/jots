#pragma once

#include <mpi.h>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/foreach.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string.hpp"
#include "mfem/mfem.hpp"
#include "precice/SolverInterface.hpp"

#include "option_structure.hpp"
#include "boundary_condition.hpp"
#include "conductivity_model.hpp"


class Config
{
    private:
        
        boost::property_tree::ptree property_tree;

        std::string input_file;            /*!< \brief Input file to parse */

        std::string mesh_file;             /*!< \brief Mesh file to read in */
        int fe_order;                 /*!< \brief FE Order (solution mapping order, not necessarily same as geometric mapping order from mesh file) */
        int serial_refine;            /*!< \brief Number of times to refine mesh before parallel decomposition */
        int parallel_refine;          /*!< \brief Number of times to refine mesh after parallel decomposition  */
        
        double density;
        double Cp;
        ConductivityModel* cond_model;                 /*!< \brief Thermal conductivity model of material */

        bool use_restart;             /*!< \brief Boolean indicating if restart file should be loaded up as initial condition */
        std::string restart_file;          /*!< \brief Restart file to load + use; only read if use_restart is true */
        double initial_temp;          /*!< \brief Initial temperature field to set; only used if use_restart is false */

        bool with_preCICE;
        std::string preCICE_participant_name;
        std::string preCICE_config_file;
        std::string preCICE_mesh_name;

        BoundaryCondition** boundary_conditions; /*!< \brief Array containing BoundaryConditions objects, value = ptr to BoundaryCondition object */
        size_t bc_count;        

        TIME_SCHEME time_scheme;      /*!< \brief Time integration scheme to use */
        double t0;                    /*!< \brief Starting time */
        double tf;                    /*!< \brief Final time to run to */
        double dt;                    /*!< \brief Delta time, timestep */

        SOLVER solver;                /*!< \brief Linear system solver type */
        PRECONDITIONER prec;          /*!< \brief Preconditioner to use */
        double abs_tol;               /*!< \brief Solver absolute tolerance */
        double rel_tol;               /*!< \brief Solver relative tolerance */
        int max_iter;                 /*!< \brief Maximum solver iterations */

        int restart_freq;             /*!< \brief Frequency to output restart files (iterations per output) */
        int vis_freq;                 /*!< \brief Frequency to output Paraview files (iterations per output) */

    public:
        Config(const char* in_file);

        // Note that these Read functions below are designed to be called in this order
        void ReadFESetup();
        void ReadAndInitMatProps();
        void ReadIC();
        void ReadpreCICE();
        void ReadAndInitBCs(mfem::ParGridFunction* in_T_gf=nullptr, precice::SolverInterface* interface=nullptr);
        void ReadTimeInt();
        void ReadLinSolSettings();
        void ReadOutput();

        std::string GetMeshFile() const { return mesh_file; }

        int GetFEOrder() const { return fe_order; }

        int GetSerialRefine() const { return serial_refine; }

        int GetParallelRefine() const { return parallel_refine; }

        double GetDensity() const { return density; }
        
        double GetCp() const { return Cp; }

        ConductivityModel* GetConductivityModel() const { return cond_model; }

        bool UsesRestart() const { return use_restart; }

        std::string GetRestartFile() const { return restart_file; }

        double GetInitialTemp() const { return initial_temp; } 

        BoundaryCondition** GetBCs() const { return boundary_conditions; }

        int GetBCCount() const {return bc_count;};

        bool UsingpreCICE() const { return with_preCICE; };

        std::string GetpreCICEParticipantName() const { return preCICE_participant_name; };
        
        std::string GetpreCICEConfigFile() const { return preCICE_config_file; };

        std::string GetpreCICEMeshName() const { return preCICE_mesh_name; };

        mfem::ODESolver* GetODESolver() const; // Returns ODESolver that must be deleted by caller!
        
        std::string GetTimeSchemeString() const;

        double GetStartTime() const { return t0; }

        double GetFinalTime() const { return tf; }

        double Getdt() const { return dt; }

        mfem::IterativeSolver* GetSolver(MPI_Comm comm_) const; // Returns IterativeSolver that must be deleted by caller!

        std::string GetSolverString() const;

        mfem::HypreSmoother::Type GetPrec() const;

        std::string GetPrecString()  const;

        int GetMaxIter() const { return max_iter; }

        double GetAbsTol() const { return abs_tol; }

        double GetRelTol() const { return rel_tol; }

        int GetRestartFreq() const { return restart_freq; }

        int GetVisFreq() const { return vis_freq; }

        void ReorderBCs(const mfem::Array<int> bdr_attributes);

        ~Config();

    protected:

};