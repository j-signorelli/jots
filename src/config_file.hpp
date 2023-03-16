#pragma once

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"

#include "option_structure.hpp"
#include "boundary_condition.hpp"
#include "conductivity_model.hpp"
class Config
{
    private:
        
        std::string input_file;            /*!< \brief Input file to parse */

        std::string mesh_file;             /*!< \brief Mesh file to read in */
        int fe_order;                 /*!< \brief FE Order (solution mapping order, not necessarily same as geometric mapping order from mesh file) */
        int serial_refine;            /*!< \brief Number of times to refine mesh before parallel decomposition */
        int parallel_refine;          /*!< \brief Number of times to refine mesh after parallel decomposition  */
        
        double density;
        double Cp;
        ConductivityModel* cond_model;                 /*!< \brief Thermal diffusivity model of material */

        bool use_restart;             /*!< \brief Boolean indicating if restart file should be loaded up as initial condition */
        std::string restart_file;          /*!< \brief Restart file to load + use; only read if use_restart is true */
        double initial_temp;          /*!< \brief Initial temperature field to set; only used if use_restart is false */

        BoundaryCondition** boundary_conditions; /*!< \brief Array where index = boundary attribute - 1, value = BoundaryCondition object */
        size_t bc_count;

        TIME_SCHEME time_scheme;      /*!< \brief Time integration scheme to use */
        double tf;                    /*!< \brief Final time to run to */
        double dt;                    /*!< \brief Delta time, timestep */

        int restart_freq;             /*!< \brief Frequency to output restart files (iterations per output) */
        int vis_freq;                 /*!< \brief Frequency to output Paraview files (iterations per output) */
        
    
    public:
        Config(const char* in_file);

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

        int GetBCCount() {return bc_count;};

        TIME_SCHEME GetTimeScheme() const { return time_scheme; }
        
        double GetFinalTime() const { return tf; }

        double Getdt() const { return dt; }

        std::string ToString() const;

        ~Config();

    protected:

};