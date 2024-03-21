
#include "config_file.hpp"

namespace bp = boost::property_tree;
using namespace std;

Config::Config(const char* in_file) : input_file(in_file)
{   
    // Parse ini file, setting property_tree
    bp::read_ini(input_file, property_tree);

    // Read in everything (setting default values if needed):
    ReadFESetup();
    ReadMatProps();
    ReadPrecice();
    ReadBCs();
    ReadAdditionalSettings();
    ReadTimeInt();
    ReadLinSolSettings();
    ReadNewtonSettings();
    ReadOutput();
}

void Config::SetInputStringVector(string in, vector<string>& output) // Comma delineated string --> no whitespace string vector
{
    boost::algorithm::split(output, in, boost::algorithm::is_any_of(","));
    
    for (size_t i = 0; i < output.size(); i++)
        boost::algorithm::trim(output[i]);// Trim whitespaces
    
}

void Config::ReadFESetup()
{   
    // Read FiniteElementSetup
    sim_type_label = property_tree.get("FiniteElementSetup.Simulation_Type", "Unsteady");
    BINARY_CHOICE restart_choice = Binary_Choice_Map.at(property_tree.get("FiniteElementSetup.Use_Restart", "No"));
    use_restart = restart_choice == BINARY_CHOICE::YES ? true : false;

    mesh_file = property_tree.get("FiniteElementSetup.Mesh_File", "mesh_name.mesh");
    fe_order = property_tree.get("FiniteElementSetup.FE_Order", 1);
    serial_refine = property_tree.get("FiniteElementSetup.Serial_Refine", 0);
    parallel_refine = property_tree.get("FiniteElementSetup.Parallel_Refine", 0);
    restart_prefix = property_tree.get("FiniteElementSetup.Restart_Prefix", "restart");
    restart_cycle = property_tree.get<int>("FiniteElementSetup.Restart_Cycle", 0);

    // Look for any solution initializations
    BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child("FiniteElementSetup"))
    {
        const string &label = v.first;
        size_t found = label.find("_Initialization");

        if (found == string::npos)
            continue;

        initialization_map[label.substr(0, found)] = boost::lexical_cast<double>(v.second.data());

        
    }
}

void Config::ReadMatProps()
{
    // Read MaterialProperties
    BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child("MaterialProperties"))
    {   
        // Get the material property label
        string label = v.first;

        // Get material property info:
        vector<string> mat_prop_info;
        SetInputStringVector(v.second.data(), mat_prop_info);
        
        // Add to mat_prop_info
        mat_prop_info_map[label] = mat_prop_info;
        
    }
}

void Config::ReadPrecice()
{
    // If no preCICE section detected, then do nothing + set with_preCICE
    if (property_tree.find("preCICE") == property_tree.not_found())
    {
        with_precice = false;
        precice_participant_name = "";
        precice_config_file = "";
    }
    else // else preCICE section exists, read in info
    {
        with_precice = true;
        precice_participant_name = property_tree.get("preCICE.Participant_Name", "jots");
        precice_config_file = property_tree.get("preCICE.Config_File", "../precice-config.xml");
    }
}

void Config::ReadBCs()
{
    // Read any section with "BoundaryConditions" existing in its header

    for (bp::ptree::iterator it = property_tree.begin(); it != property_tree.end(); it++)
    {
        const string &root = it->first;
        size_t found = root.find("BoundaryConditions");

        if (found == string::npos)
            continue;
        
        string type = root.substr(0, found);

        stringstream label;
        label << type << "BoundaryConditions";
        // Read "[__type__BoundaryConditions]"
        BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child(label.str()))
        {   
            // Get the attribute, split by "_" and get the final element
            vector<string> label;
            boost::algorithm::split(label, v.first, boost::algorithm::is_any_of("_"));

            string s_attr = label.back();
            int attr = stoi(s_attr);

            // Get the BC type:
            vector<string> single_bc_info;
            SetInputStringVector(v.second.data(), single_bc_info);
            
            // Add to bc_info_map
            bc_info_map[type][attr] = single_bc_info;
            
        }
    }

}

void Config::ReadAdditionalSettings()
{
    if (property_tree.find("AdditionalSettings") != property_tree.not_found())
    {
        // Read "[AdditionalSettings]"
        BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child("AdditionalSettings"))
        {   
            // Get the setting and value
            string setting = v.first.data();
            string value = v.second.data();

            // Add to map
            additional_settings[setting] = value;
                
        }
    }
}

void Config::ReadTimeInt()
{
    // Read TimeIntegration
    // If no time integration, then do nothing and set using_time_integration
    if (property_tree.find("TimeIntegration") == property_tree.not_found())
    {
        using_time_integration = false;
        time_scheme_label = "Euler_Implicit";
        dt = 0.1;
        max_timesteps = 100;
        time_print_freq = 1;
    }
    else
    {
        using_time_integration = true;
        time_scheme_label = property_tree.get("TimeIntegration.Time_Scheme", "Euler_Implicit");
        dt = property_tree.get("TimeIntegration.Delta_Time", 0.1);
        max_timesteps = property_tree.get("TimeIntegration.Max_Timesteps", 100);
        time_print_freq = property_tree.get("TimeIntegration.Print_Freq", 1);
    }
}

void Config::ReadLinSolSettings()
{
    // Read LinearSolverSettings
    solver_label = property_tree.get("LinearSolverSettings.Solver", "FGMRES");
    prec_label = property_tree.get("LinearSolverSettings.Preconditioner", "Chebyshev");
    abs_tol = property_tree.get("LinearSolverSettings.Absolute_Tolerance", 1e-16);
    rel_tol = property_tree.get("LinearSolverSettings.Relative_Tolerance", 1e-10);
    max_iter = property_tree.get("LinearSolverSettings.Max_Iterations", 100);
    string print_level = property_tree.get("LinearSolverSettings.Print_Level", "Errors, Warnings");
    SetInputStringVector(print_level, ls_print_level);
}

void Config::ReadNewtonSettings()
{
    // Read NewtonSolverSettings
    if (property_tree.find("NewtonSolverSettings") == property_tree.not_found())
    {
        using_newton = false;
        newton_max_iter = 0;
        newton_abs_tol = 1e-16;
        newton_rel_tol = 1e-10;
        SetInputStringVector("Errors, Warnings, Iterations", newton_print_level);
    }
    else
    {
        using_newton=true;
        newton_max_iter = property_tree.get("NewtonSolverSettings.Max_Iterations", 100);
        newton_abs_tol = property_tree.get("NewtonSolverSettings.Absolute_Tolerance", 1e-16);
        newton_rel_tol = property_tree.get("NewtonSolverSettings.Relative_Tolerance", 1e-10);
        string print_level = property_tree.get("NewtonSolverSettings.Print_Level", "Errors, Warnings, Iterations");
        SetInputStringVector(print_level, newton_print_level);
    }
}

void Config::ReadOutput()
{
    //Read Output
    restart_freq = property_tree.get("Output.Restart_Freq", 10);
    vis_freq = property_tree.get("Output.Visualization_Freq", 10);

}