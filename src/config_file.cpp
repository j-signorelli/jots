#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"

#include "config_file.hpp"

namespace bp = boost::property_tree;
using namespace std;

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

    // Read InitialCondition
    BINARY_CHOICE restart_choice = Binary_Choice_Map.at(property_tree.get<string>("InitialCondition.Use_Restart"));
    m_use_restart = restart_choice == BINARY_CHOICE::YES ? true : false;

    if (!m_use_restart)
    {
        // Grab initial temperature field set
        m_initial_temp = property_tree.get<double>("InitialCondition.Initial_Temperature");
    }
    else
    {
        // Else need to read in restart file
    }
    // Read BoundaryConditions
    int bc_count = property_tree.get_child("BoundaryConditions").size()/2;
    for (int i = 0; i < bc_count; i++)
    {
        // Insert type and value
        BOUNDARY_CONDITION type = Boundary_Condition_Map.at(property_tree.get<string>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Type"));
        double value =  property_tree.get<double>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Value");
        m_boundary_conditions[i+1] = make_tuple(type, value);
    }

    // Read TimeIntegration
    m_time_scheme = Time_Scheme_Map.at(property_tree.get<string>("TimeIntegration.Time_Scheme"));
    m_dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    m_tf = property_tree.get<double>("TimeIntegration.Final_Time");

    //Read Output
    m_restart_freq = property_tree.get<int>("Output.Restart_Freq");
    m_vis_freq = property_tree.get<int>("Output.Visualization_Freq");


}

string Config::ToString() const
{
    string s = "";
    s += "Mesh File: " + m_mesh_file + "\n";
    s += "FE Order: " + to_string(m_fe_order) + "\n";
    s += "Serial Refine: " + to_string(m_serial_refine) + "\n";
    s += "Parallel Refine: " + to_string(m_parallel_refine) + "\n";
    s += "Kappa: " + to_string(m_kappa) + "\n";
    
    for (int i = 0; i < m_boundary_conditions.size(); i++)
    {
        int type = static_cast<int>(get<0>(m_boundary_conditions.at(i+1))); //.at is used here instead of [] as [] is non const, while this fxn is meant to be const
        double value = get<1>(m_boundary_conditions.at(i+1));
        s += "Boundary " + to_string(i+1) + " Type, Value: " + to_string(type) + ", " + to_string(value) + "\n";
    }
    s += "Time Scheme: " + to_string(static_cast<int>(m_time_scheme)) + "\n";
    s += "dt: " + to_string(m_dt) + "\n";
    s += "tf: " + to_string(m_tf) + "\n";
    s += "Restart Freq: " + to_string(m_restart_freq) + "\n";
    s += "Visualization Freq: " + to_string(m_vis_freq) + "\n";

    return s;
}