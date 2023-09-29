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

        SIMULATION_TYPE sim_type;
        std::string mesh_file;             /*!< \brief Mesh file to read in */
        int fe_order;                 /*!< \brief FE Order (solution mapping order, not necessarily same as geometric mapping order from mesh file) */
        int serial_refine;            /*!< \brief Number of times to refine mesh before parallel decomposition */
        int parallel_refine;          /*!< \brief Number of times to refine mesh after parallel decomposition  */
        
        // Material properties:
        double density;
        std::vector<std::string> specific_heat_info;
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

        std::string GetInputFile() const { return input_file; };

        void SetInputFile(std::string in_file) { input_file = in_file; };

        SIMULATION_TYPE GetSimType() const { return sim_type; };

        void SetSimType(SIMULATION in_type) { sim_type = in_type; };

        std::string GetMeshFile() const { return mesh_file; };

        void SetMeshFile(std::string in_file) { mesh_file = in_file; };

        int GetFEOrder() const { return fe_order; };

        void SetFEOrder(int in_order) { fe_order = in_order; };
        
        int GetSerialRefine() const { return serial_refine; };

        void SetSerialRefine(int in_ref) {serial_refine = in_ref; };

        int GetParallelRefine() const { return parallel_refine; };

        void SetParallelRefine(int in_ref) { parallel_refine = in_ref; };

        double GetDensity() const { return density; };
        
        void SetDensity(double in_rho) { density = in_rho; };

        std::vector<std::string> GetSpecificHeatInfo() const { return specific_heat_info; };

        void SetSpecificHeatInfo(std::vector<std::string> in_info) { specific_heat_info = in_info; };

        std::vector<std::string> GetCondInfo() const { return conductivity_info; };
        
        void SetCondInfo(std::vector<std::string> in_info) { conductivity_info = in_info; };

        bool UsesRestart() const { return use_restart; };

        void SetRestart(bool in_restart) { use_restart = in_restart; };

        std::string GetRestartPrefix() const { return restart_prefix; };
        
        void SetRestartPrefix(std::string in_prefix) { restart_prefix = in_prefix; };

        int GetRestartCycle() const { return restart_cycle; };

        void SetRestartCycle(int in_cycle) { restart_cycle = in_cycle; };

        double GetInitialTemp() const { return initial_temp; } ;

        void SetInitialTemp(double in_temp) { initial_temp = in_temp; };

        int GetBCCount() const {return bc_info.size(); };

        void ResizeBCs(int in_bc_count) { bc_info.resize(in_bc_count); };

        bool UsingPrecice() const { return with_precice; };

        void SetPrecice(bool in_precice) { with_precice = in_precice; };

        std::string GetPreciceParticipantName() const { return precice_participant_name; };
        
        void SetPreciceParticipantName(std::string in_name) { precice_participant_name = in_name; };

        std::string GetPreciceConfigFile() const { return precice_config_file; };
    
        void SetPreciceConfigFile(std::string in_file) { precice_config_file = in_file; };

        std::pair<int, std::vector<std::string>> GetBCInfo(int i) const { return bc_info[i]; };

        void SetBCInfo(int i, std::pair<int, std::vector<std::string>> in_bc) { bc_info[i] = in_bc; };

        mfem::ODESolver* GetODESolver() const; // Returns ODESolver that must be deleted by caller!
        
        void SetODESolver(TIME_SCHEME in_scheme) { time_scheme = in_scheme; };

        std::string GetTimeSchemeString() const;

        double GetFinalTime() const { return tf; };

        void SetFinalTime(double in_tf) { tf = in_tf; };

        double Getdt() const { return dt; };

        void Setdt(double in_dt) { dt = in_dt; };

        mfem::IterativeSolver* GetSolver(MPI_Comm comm_) const; // Returns IterativeSolver that must be deleted by caller!

        void SetSolver(SOLVER in_solver) { solver = in_solver; };

        std::string GetSolverString() const;

        mfem::HypreSmoother::Type GetPrec() const;

        void SetPrec(PRECONDITIONER in_prec) { prec = in_prec; };

        std::string GetPrecString()  const;

        int GetMaxIter() const { return max_iter; };

        void SetMaxIter(int in_iter) { max_iter = in_iter; };

        double GetAbsTol() const { return abs_tol; };

        void SetAbsTol(double in_tol) { abs_tol = in_tol; };

        double GetRelTol() const { return rel_tol; };

        void SetRelTol(double in_tol) { rel_tol = in_tol; };

        int GetRestartFreq() const { return restart_freq; };

        void SetRestartFreq(int in_freq) { restart_freq = in_freq; };

        int GetVisFreq() const { return vis_freq; };

        void SetVisFreq(int in_freq) { vis_freq = in_freq; };

    protected:

};