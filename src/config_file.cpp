
#include "config_file.hpp"

namespace bp = boost::property_tree;
using namespace std;

Config::Config(const char* in_file) : input_file(in_file)
{   

    // Set up property_tree
    bp::ptree property_tree;

    // Parse ini file
    bp::read_ini(input_file, property_tree);

    // Set all member variables using values from tree

    // Read FiniteElementSetup
    mesh_file = property_tree.get<string>("FiniteElementSetup.Mesh_File");
    fe_order = property_tree.get<int>("FiniteElementSetup.FE_Order");
    serial_refine = property_tree.get<int>("FiniteElementSetup.Serial_Refine");
    parallel_refine = property_tree.get<int>("FiniteElementSetup.Parallel_Refine");


    // Read MaterialProperties
    density = property_tree.get<double>("MaterialProperties.Density");
    Cp = property_tree.get<double>("MaterialProperties.Specific_Heat_Cp");
    CONDUCTIVITY_MODEL model = Conductivity_Model_Map.at(property_tree.get<string>("MaterialProperties.Thermal_Conductivity_Model"));
    switch (model)
    {
        case CONDUCTIVITY_MODEL::CONSTANT:
            //double kappa = property_tree.get<double>("MaterialProperties.Kappa");
            cond_model = new ConstantCond(property_tree.get<double>("MaterialProperties.k"));
            break;
        case CONDUCTIVITY_MODEL::LINEARIZED:
            cond_model = new LinearizedCond(property_tree.get<double>("MaterialProperties.k"), property_tree.get<double>("MaterialProperties.alpha"));
            break;
    }

    // Read InitialCondition
    BINARY_CHOICE restart_choice = Binary_Choice_Map.at(property_tree.get<string>("InitialCondition.Use_Restart"));
    use_restart = restart_choice == BINARY_CHOICE::YES ? true : false;

    if (!use_restart)
    {
        // Grab initial temperature field set
        initial_temp = property_tree.get<double>("InitialCondition.Initial_Temperature");
    }
    else
    {
        // Else need to read in restart file - save file name
        restart_file = property_tree.get<string>("InitialCondition.Restart_File");
    }
    // Read BoundaryConditions
    int bc_count = property_tree.get_child("BoundaryConditions").size()/2;

    for (int i = 0; i < bc_count; i++)
    {
        // For now: just assume each BC type has a type and a value. May update in future
        BOUNDARY_CONDITION type = Boundary_Condition_Map.at(property_tree.get<string>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Type"));
        double value =  property_tree.get<double>("BoundaryConditions.Boundary_" + to_string(i+1)  + "_Value");
        
        switch (type)
        {
            case BOUNDARY_CONDITION::HEATFLUX:
                boundary_conditions.push_back(new UniformHeatFluxBC(value));
                break;
            case BOUNDARY_CONDITION::ISOTHERMAL:
                boundary_conditions.push_back(new UniformIsothermalBC(value));
                break;

        }
        
    }

    // Read TimeIntegration
    time_scheme = Time_Scheme_Map.at(property_tree.get<string>("TimeIntegration.Time_Scheme"));
    dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    tf = property_tree.get<double>("TimeIntegration.Final_Time");

    //Read Output
    restart_freq = property_tree.get<int>("Output.Restart_Freq");
    vis_freq = property_tree.get<int>("Output.Visualization_Freq");


}

/* OLD: used to debug
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
*/