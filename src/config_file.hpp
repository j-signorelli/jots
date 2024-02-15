#pragma once

#include <mpi.h>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/foreach.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string.hpp"

#include "option_structure.hpp"
#include "jots_common.hpp"

class Config
{
    private:
        
        boost::property_tree::ptree property_tree;

        std::string input_file;            /*!< \brief Input file to parse */

        std::string sim_type_label;
        std::string mesh_file;             /*!< \brief Mesh file to read in */
        int fe_order;                 /*!< \brief FE Order (solution mapping order, not necessarily same as geometric mapping order from mesh file) */
        int serial_refine;            /*!< \brief Number of times to refine mesh before parallel decomposition */
        int parallel_refine;          /*!< \brief Number of times to refine mesh after parallel decomposition  */
        
        // Material properties:
        std::map<std::string, std::vector<std::string>> mat_prop_info_map; // Map of material properties, key is mat prop label and value is string vector for that mat prop

        bool use_restart;             /*!< \brief Boolean indicating if restart file should be loaded up as initial condition */
        std::string restart_prefix;          /*!< \brief Restart file to load + use; only read if use_restart is true */
        int restart_cycle;
        double initial_temp;          /*!< \brief Initial temperature field to set; only used if use_restart is false */

        bool with_precice;
        std::string precice_participant_name;
        std::string precice_config_file;

        size_t bc_count;     
        std::map<int, std::vector<std::string>> bc_info_map; // Map where key is attribute, value is string vector for that BC
        
        bool using_time_integration;
        std::string time_scheme_label;      /*!< \brief Time integration scheme to use */
        int max_timesteps;             /*!< \brief Delta time, timestep */
        double dt;                    /*!< \brief Delta time, timestep */
        int time_print_freq;

        std::string solver_label;                /*!< \brief Linear system solver type */
        std::string prec_label;          /*!< \brief Preconditioner to use */
        double abs_tol;               /*!< \brief Solver absolute tolerance */
        double rel_tol;               /*!< \brief Solver relative tolerance */
        int max_iter;                 /*!< \brief Maximum solver iterations */
        
        bool using_newton;
        double newton_abs_tol;               /*!< \brief Solver absolute tolerance */
        double newton_rel_tol;               /*!< \brief Solver relative tolerance */
        int newton_max_iter;                 /*!< \brief Maximum solver iterations */

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
        void ReadNewtonSettings();
        void ReadOutput();

    public:
        Config(const char* in_file);

        Config(const std::string in_file) : Config(in_file.c_str()) {};

        std::string GetInputFile() const { return input_file; };

        void SetInputFile(std::string in_file) { input_file = in_file; };

        std::string GetSimTypeLabel() const { return sim_type_label; };

        void SetSimTypeLabel(std::string in_type_label) { sim_type_label = in_type_label; };

        std::string GetMeshFile() const { return mesh_file; };

        void SetMeshFile(std::string in_file) { mesh_file = in_file; };

        int GetFEOrder() const { return fe_order; };

        void SetFEOrder(int in_order) { fe_order = in_order; };
        
        int GetSerialRefine() const { return serial_refine; };

        void SetSerialRefine(int in_ref) {serial_refine = in_ref; };

        int GetParallelRefine() const { return parallel_refine; };

        void SetParallelRefine(int in_ref) { parallel_refine = in_ref; };

        std::vector<std::string> GetMaterialPropertyKeys() const { return Helper::GetKeyVector(mat_prop_info_map); };

        std::vector<std::string> GetMaterialPropertyInfo(std::string mat_prop_label) const { return mat_prop_info_map.at(mat_prop_label); };

        void SetMaterialPropertyInfo(std::string mat_prop_label, std::vector<std::string> info) { mat_prop_info_map[mat_prop_label] = info; };

        bool UsesRestart() const { return use_restart; };

        void SetRestart(bool in_restart) { use_restart = in_restart; };

        std::string GetRestartPrefix() const { return restart_prefix; };
        
        void SetRestartPrefix(std::string in_prefix) { restart_prefix = in_prefix; };

        int GetRestartCycle() const { return restart_cycle; };

        void SetRestartCycle(int in_cycle) { restart_cycle = in_cycle; };

        double GetInitialTemp() const { return initial_temp; } ;

        void SetInitialTemp(double in_temp) { initial_temp = in_temp; };

        int GetBCCount() const {return bc_info_map.size(); };

        bool UsingPrecice() const { return with_precice; };

        void SetPrecice(bool in_precice) { with_precice = in_precice; };

        std::string GetPreciceParticipantName() const { return precice_participant_name; };
        
        void SetPreciceParticipantName(std::string in_name) { precice_participant_name = in_name; };

        std::string GetPreciceConfigFile() const { return precice_config_file; };
    
        void SetPreciceConfigFile(std::string in_file) { precice_config_file = in_file; };

        std::vector<int> GetBCKeys() const { return Helper::GetKeyVector(bc_info_map); };

        std::vector<std::string> GetBCInfo(int attr) const { return bc_info_map.at(attr); };

        void SetBCInfo(int attr, std::vector<std::string> in_bc) { bc_info_map[attr] = in_bc; };

        void DeleteBCInfo(int attr) { bc_info_map.erase(attr); };
        
        bool UsingTimeIntegration() const { return using_time_integration; };

        void SetTimeIntegration(bool in_using) { using_time_integration = in_using; };

        void SetTimeSchemeLabel(std::string in_scheme) { time_scheme_label = in_scheme; };

        std::string GetTimeSchemeLabel() const { return time_scheme_label; };

        int GetMaxTimesteps() const { return max_timesteps; };

        void SetMaxTimesteps(int in_max_timesteps) {  max_timesteps = in_max_timesteps; };

        double Getdt() const { return dt; };

        void Setdt(double in_dt) { dt = in_dt; };

        int GetTimePrintFreq() const { return time_print_freq; };

        void SetTimePrintFreq(int in_freq) { time_print_freq = in_freq; };
        
        std::string GetSolverLabel() const { return solver_label; };

        void SetSolverLabel(std::string in_solver) { solver_label = in_solver; };

        std::string GetPrecLabel() const { return prec_label; };

        void SetPrecLabel(std::string in_prec) { prec_label = in_prec; };

        int GetMaxIter() const { return max_iter; };

        void SetMaxIter(int in_iter) { max_iter = in_iter; };

        double GetAbsTol() const { return abs_tol; };

        void SetAbsTol(double in_tol) { abs_tol = in_tol; };

        double GetRelTol() const { return rel_tol; };

        void SetRelTol(double in_tol) { rel_tol = in_tol; };

        bool UsingNewton() const { return using_newton; };
        
        void SetNewton(bool in_newton) { using_newton = in_newton; };

        int GetNewtonMaxIter() const { return newton_max_iter; };

        void SetNewtonMaxIter(int in_iter) { newton_max_iter = in_iter; };

        double GetNewtonAbsTol() const { return newton_abs_tol; };

        void SetNewtonAbsTol(double in_tol) { newton_abs_tol = in_tol; };

        double GetNewtonRelTol() const { return newton_rel_tol; };

        void SetNewtonRelTol(double in_tol) { newton_rel_tol = in_tol; };

        int GetRestartFreq() const { return restart_freq; };

        void SetRestartFreq(int in_freq) { restart_freq = in_freq; };

        int GetVisFreq() const { return vis_freq; };

        void SetVisFreq(int in_freq) { vis_freq = in_freq; };

    protected:

};