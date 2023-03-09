#include "config_file.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include <iostream>

#include "boost/filesystem.hpp"

using namespace std;

namespace bp = boost::property_tree;

Config::Config(const char* input_file) : m_input_file(input_file)
{
    // Set up property_tree
    bp::ptree property_tree;

    // Parse ini file
    bp::read_ini(m_input_file, property_tree);

    // Set all member variables using values from tree
    m_mesh_file = property_tree.get<string>("FiniteElementSetup.Mesh_File");
    m_fe_order = property_tree.get<int>("FiniteElementSetup.FE_Order");
    m_serial_refine = property_tree.get<int>("FiniteElementSetup.Serial_Refine");
    m_parallel_refine = property_tree.get<int>("FiniteElementSetup.Parallel_Refine");

    m_kappa = property_tree.get<double>("MaterialProperties.Thermal_Diffusivity");

    m_time_scheme = property_tree.get<int>("TimeIntegration.Time_Scheme");
    m_dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    m_tf = property_tree.get<double>("TimeIntegration.Final_Time");

    m_restart_freq = property_tree.get<int>("Output.Restart_Freq");
    m_vis_freq = property_tree.get<int>("Output.Visualization_Freq");


}
