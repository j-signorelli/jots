
#include "config_file.hpp"

namespace bp = boost::property_tree;
using namespace std;
using namespace mfem;
using namespace precice;

Config::Config(const char* in_file) : input_file(in_file)
{   
    // Parse ini file, setting property_tree
    bp::read_ini(input_file, property_tree);

    // Read in everything (setting default values if needed):
    ReadFESetup();
    ReadMatProps();
    ReadPrecice();
    ReadBCs();
    ReadTimeInt();
    ReadLinSolSettings();
    ReadOutput();
}

void Config::SetInputStringVector(string in, vector<string>& output) // Comma delineated string --> no whitespace string vector
{
    boost::algorithm::split(output, in, boost::algorithm::is_any_of(","));
    
    for (int i = 0; i < output.size(); i++)
        boost::algorithm::trim(output[i]);// Trim whitespaces
    
}

void Config::ReadFESetup()
{   
    // Read FiniteElementSetup
    sim_type = Simulation_Type_Map.at(property_tree.get("FiniteElementSetup.Simulation_Type", SIMULATION_TYPE::UNSTEADY))
    BINARY_CHOICE restart_choice = Binary_Choice_Map.at(property_tree.get("FiniteElementSetup.Use_Restart", "No"));
    use_restart = restart_choice == BINARY_CHOICE::YES ? true : false;

    mesh_file = property_tree.get("FiniteElementSetup.Mesh_File", "mesh_name.mesh");
    fe_order = property_tree.get("FiniteElementSetup.FE_Order", 1);
    serial_refine = property_tree.get("FiniteElementSetup.Serial_Refine", 0);
    parallel_refine = property_tree.get("FiniteElementSetup.Parallel_Refine", 0);
    initial_temp = property_tree.get<double>("FiniteElementSetup.Initial_Temperature", 100.0);
    restart_prefix = property_tree.get("FiniteElementSetup.Restart_Prefix", "restart");
    restart_cycle = property_tree.get<int>("FiniteElementSetup.Restart_Cycle", 0);

}

void Config::ReadMatProps()
{
    // Read MaterialProperties
    density = property_tree.get("MaterialProperties.Density", 1.0);
    SetInputStringVector(property_tree.get("MaterialProperties.Specific_Heat_C", "Uniform, 1000"), specific_heat_info);
    SetInputStringVector(property_tree.get("MaterialProperties.Thermal_Conductivity_k", "Uniform, 100"), conductivity_info);
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
    // Read BoundaryConditions
    bc_count = property_tree.get_child("BoundaryConditions").size();

    BOOST_FOREACH(const bp::ptree::value_type &v , property_tree.get_child("BoundaryConditions"))
    {   
        // Get the attribute, split by "_" and get the final element
        vector<string> label;
        boost::algorithm::split(label, v.first, boost::algorithm::is_any_of("_"));

        string s_attr = label.back();
        int attr = stoi(s_attr);

        // Get the BC type:
        vector<string> single_bc_info;
        SetInputStringVector(v.second.data(), single_bc_info);
        
        // Add to bc_info
        bc_info.push_back(make_pair(attr, single_bc_info));
        
    }
}

void Config::ReadTimeInt()
{
    // Read TimeIntegration
    time_scheme = Time_Scheme_Map.at(property_tree.get("TimeIntegration.Time_Scheme", "Euler_Implicit"));
    dt = property_tree.get("TimeIntegration.Delta_Time", 0.1);
    tf = property_tree.get("TimeIntegration.Final_Time", 10.0);
}

void Config::ReadLinSolSettings()
{
    // Read LinearSolverSettings
    solver = Solver_Map.at(property_tree.get("LinearSolverSettings.Solver", "FGMRES"));
    prec = Preconditioner_Map.at(property_tree.get("LinearSolverSettings.Preconditioner", "Chebyshev"));
    abs_tol = property_tree.get("LinearSolverSettings.Absolute_Tolerance", 1e-16);
    rel_tol = property_tree.get("LinearSolverSettings.Relative_Tolerance", 1e-10);
    max_iter = property_tree.get("LinearSolverSettings.Max_Iterations", 100);
}

void Config::ReadOutput()
{
    //Read Output
    restart_freq = property_tree.get("Output.Restart_Freq", 10);
    vis_freq = property_tree.get("Output.Visualization_Freq", 10);

}


ODESolver* Config::GetODESolver() const
{
    switch (time_scheme)
    {
        case TIME_SCHEME::EULER_IMPLICIT:
            return new BackwardEulerSolver;
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            return new ForwardEulerSolver;
            break;
        case TIME_SCHEME::RK4:
            return new RK4Solver;
    }
}

string Config::GetTimeSchemeString() const
{
    switch (time_scheme)
    {
        case TIME_SCHEME::EULER_IMPLICIT:
            return "Euler Implicit";
            break;
        case TIME_SCHEME::EULER_EXPLICIT:
            return "Euler Explicit";
            break;
        case TIME_SCHEME::RK4:
            return "RK4";
    }

}

IterativeSolver* Config::GetSolver(MPI_Comm comm_) const
{
    switch (solver)
    {
        case SOLVER::CG:
            return new CGSolver(comm_);
            break;
        case SOLVER::GMRES:
            return new GMRESSolver(comm_);
            break;
        case SOLVER::FGMRES:
            return new FGMRESSolver(comm_);
            break;
    }
}

string Config::GetSolverString() const
{
    switch (solver)
    {
        case SOLVER::CG:
            return "Conjugate Gradient";
            break;
        case SOLVER::GMRES:
            return "GMRES";
            break;
        case SOLVER::FGMRES:
            return "FGMRES";
            break;
    }

}
// TODO: clean up these "GetStrings" such that not needing to copy + paste ideally
HypreSmoother::Type Config::GetPrec() const
{
    switch (prec)
    {
        case PRECONDITIONER::JACOBI:
            return HypreSmoother::Jacobi;
            break;
        case PRECONDITIONER::CHEBYSHEV:
            return HypreSmoother::Chebyshev;
            break;
    }
}

string Config::GetPrecString() const
{
    switch (prec)
    {
        case PRECONDITIONER::JACOBI:
            return "Jacobi";
            break;
        case PRECONDITIONER::CHEBYSHEV:
            return "Chebyshev";
            break;
    }
}
