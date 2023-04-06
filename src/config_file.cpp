
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
    bc_count = property_tree.get_child("BoundaryConditions").size();
    boundary_conditions = new BoundaryCondition*[bc_count];
    int index = 0;
    BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child("BoundaryConditions"))
    {   
        // Get the attribute:
        char s_attr = v.first.back();
        int attr = atoi(&s_attr);

        // Get the BC type:
        vector<string> bc_info;
        boost::algorithm::split(bc_info, v.second.data(), boost::algorithm::is_any_of(","));
        
        // Clear out any leading or trailing white space
        boost::algorithm::trim(bc_info[0]);
        boost::algorithm::trim(bc_info[1]);

        // Get the value set
        double value =  stod(bc_info[1].c_str());

        // Set the BC
        BOUNDARY_CONDITION type = Boundary_Condition_Map.at(bc_info[0]);
        
        switch (type)
        {
            case BOUNDARY_CONDITION::HEATFLUX:
                boundary_conditions[index] = new UniformHeatFluxBC(attr, value);
                break;
            case BOUNDARY_CONDITION::ISOTHERMAL:
                boundary_conditions[index] = new UniformIsothermalBC(attr, value);
                break;

        }

        index++;
        
        
        
    }

    // TODO: Read LinearSystemSettings

    // Read TimeIntegration
    time_scheme = Time_Scheme_Map.at(property_tree.get<string>("TimeIntegration.Time_Scheme"));
    dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    tf = property_tree.get<double>("TimeIntegration.Final_Time");

    //Read Output
    restart_freq = property_tree.get<int>("Output.Restart_Freq");
    vis_freq = property_tree.get<int>("Output.Visualization_Freq");

    // TODO: update to automatically retrieve from restarts, but for now:
    t0 = 0;
}

Config::~Config()
{
    delete cond_model;
    for (size_t i = 0; i < bc_count; i++)
        delete boundary_conditions[i];

    delete[] boundary_conditions;
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