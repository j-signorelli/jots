#include "config_file.hpp"

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include <iostream>

using namespace std;

namespace bp = boost::property_tree;

Config::Config(const char* input_file) : m_input_file(input_file)
{
    // Set up property_tree
    bp::ptree property_tree;

    // Parse ini file
    bp::read_ini(m_input_file, property_tree);

    // Set all member variables using values from tree

    // Read FiniteElementSetup
    m_mesh_file = property_tree.get<string>("FiniteElementSetup.Mesh_File");
    m_fe_order = property_tree.get<int>("FiniteElementSetup.FE_Order");
    m_serial_refine = property_tree.get<int>("FiniteElementSetup.Serial_Refine");
    m_parallel_refine = property_tree.get<int>("FiniteElementSetup.Parallel_Refine");

    // Read MaterialProperties
    m_kappa = property_tree.get<double>("MaterialProperties.Thermal_Diffusivity");

    // Read InitialConditions
    // TODO Later -- first testing

    // Read BoundaryConditions
    // NOTE: Need to think about preCICE BC
    //          What BC to start out with? (depends on coupling scheme)
    //          Isothermal or not???
    int bc_count = property_tree.count("BoundaryConditions") / 2;

    for (int i = 0; i < bc_count; i++)
    {
        // Insert type and value
        const BOUNDARY_CONDITION type = Boundary_Condition_Map.at(property_tree.get<string>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Type"));
        const double value =  property_tree.get<double>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Value");
        m_boundary_conditions.insert(make_pair(i+1, tuple<BOUNDARY_CONDITION, double>(type, value)));
    }

    // Read TimeIntegration
    m_time_scheme = Time_Scheme_Map.at(property_tree.get<string>("TimeIntegration.Time_Scheme"));
    m_dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    m_tf = property_tree.get<double>("TimeIntegration.Final_Time");

    //Read Output
    m_restart_freq = property_tree.get<int>("Output.Restart_Freq");
    m_vis_freq = property_tree.get<int>("Output.Visualization_Freq");


}
