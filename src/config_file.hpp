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
        std::vector<std::string> conductivity_info;

        bool use_restart;             /*!< \brief Boolean indicating if restart file should be loaded up as initial condition */
        std::string restart_prefix;          /*!< \brief Restart file to load + use; only read if use_restart is true */
        int restart_cycle;
        double initial_temp;          /*!< \brief Initial temperature field to set; only used if use_restart is false */

        bool with_precice;
        std::string precice_participant_name;
        std::string precice_config_file;

        size_t bc_count;        
        std::vector<std::pair<int, std::vector<std::string>>> bc_info; // Array of pairs where first value is attribute, second is string vector for that BC
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

        void SetInputStringVector(std::string in, std::vector<std::string>& output); // Comma delineated string --> no whitespace string vector

        void ReadFESetup();
        void ReadMatProps();
        void ReadIC();
        void ReadPrecice();
        void ReadBCs();
        void ReadTimeInt();
        void ReadLinSolSettings();
        void ReadOutput();

    public:
        Config(const char* in_file);

        std::string GetMeshFile() const { return mesh_file; };

        int GetFEOrder() const { return fe_order; };

        int GetSerialRefine() const { return serial_refine; };

        int GetParallelRefine() const { return parallel_refine; };

        double GetDensity() const { return density; };
        
        double GetCp() const { return Cp; };

        std::vector<std::string> GetCondInfo() const { return conductivity_info; };
        
        bool UsesRestart() const { return use_restart; };

        std::string GetRestartPrefix() const { return restart_prefix; };
        
        int GetRestartCycle() const { return restart_cycle; };

        double GetInitialTemp() const { return initial_temp; } ;

        int GetBCCount() const {return bc_count;};

        bool UsingPrecice() const { return with_precice; };

        std::string GetPreciceParticipantName() const { return precice_participant_name; };
        
        std::string GetPreciceConfigFile() const { return precice_config_file; };
    
        std::pair<int, std::vector<std::string>> GetBCInfo(int i) const { return bc_info[i]; };

        mfem::ODESolver* GetODESolver() const; // Returns ODESolver that must be deleted by caller!
        
        std::string GetTimeSchemeString() const;

        double GetFinalTime() const { return tf; };

        double Getdt() const { return dt; };

        mfem::IterativeSolver* GetSolver(MPI_Comm comm_) const; // Returns IterativeSolver that must be deleted by caller!

        std::string GetSolverString() const;

        mfem::HypreSmoother::Type GetPrec() const;

        std::string GetPrecString()  const;

        int GetMaxIter() const { return max_iter; };

        double GetAbsTol() const { return abs_tol; };

        double GetRelTol() const { return rel_tol; };

        int GetRestartFreq() const { return restart_freq; };

        int GetVisFreq() const { return vis_freq; };

    protected:

};