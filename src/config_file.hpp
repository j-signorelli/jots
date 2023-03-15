#pragma once
#include "option_structure.hpp"

class Config
{
    private:
        
        std::string m_input_file;            /*!< \brief Input file to parse */

        std::string m_mesh_file;             /*!< \brief Mesh file to read in */
        int m_fe_order;                 /*!< \brief FE Order (solution mapping order, not necessarily same as geometric mapping order from mesh file) */
        int m_serial_refine;            /*!< \brief Number of times to refine mesh before parallel decomposition */
        int m_parallel_refine;          /*!< \brief Number of times to refine mesh after parallel decomposition  */
        
        double m_kappa;                 /*!< \brief Thermal diffusivity of material */

        bool m_use_restart;             /*!< \brief Boolean indicating if restart file should be loaded up as initial condition */
        std::string m_restart_file;          /*!< \brief Restart file to load + use; only read if m_use_restart is true */
        double m_initial_temp;          /*!< \brief Initial temperature field to set; only used if m_use_restart is false */

        std::map<int, std::tuple<BOUNDARY_CONDITION, double>> m_boundary_conditions; /*!< \brief Map where key = boundary attribute, value = (BC type, value) */

        TIME_SCHEME m_time_scheme;      /*!< \brief Time integration scheme to use */
        double m_tf;                    /*!< \brief Final time to run to */
        double m_dt;                    /*!< \brief Delta time, timestep */

        int m_restart_freq;             /*!< \brief Frequency to output restart files (iterations per output) */
        int m_vis_freq;                 /*!< \brief Frequency to output Paraview files (iterations per output) */
        
    
    public:
        Config(const char* input_file);

        std::string GetMeshFile() const 
        {
            return m_mesh_file;
        }

        int GetFEOrder() const
        {
            return m_fe_order;
        }

        int GetSerialRefine() const
        {
            return m_serial_refine;
        }

        int GetParallelRefine() const
        {
            return m_parallel_refine;
        }

        double GetKappa() const
        {
            return m_kappa;
        }

        bool UsesRestart() const
        {
            return m_use_restart;
        }

        double GetInitialTemp() const
        {
            return m_initial_temp;
        } 

        std::map<int, std::tuple<BOUNDARY_CONDITION, double>> GetBCs() const
        {
            return m_boundary_conditions;
        }

        TIME_SCHEME GetTimeScheme() const
        {
            return m_time_scheme;
        }
        
        double GetFinalTime() const
        {
            return m_tf;
        }

        double Getdt() const
        {
            return m_dt;
        }

        std::string ToString() const;


    protected:

};