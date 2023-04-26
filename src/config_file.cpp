
#include "config_file.hpp"

namespace bp = boost::property_tree;
using namespace std;
using namespace mfem;
using namespace precice;

Config::Config(const char* in_file) : input_file(in_file)
{   

    // Parse ini file
    bp::read_ini(input_file, property_tree);

    // Set all member variables using values from tree (implemented in member fxns - must be called separately)
}

void Config::ReadFESetup()
{
    // Read FiniteElementSetup
    mesh_file = property_tree.get<string>("FiniteElementSetup.Mesh_File");
    fe_order = property_tree.get<int>("FiniteElementSetup.FE_Order");
    serial_refine = property_tree.get<int>("FiniteElementSetup.Serial_Refine");
    parallel_refine = property_tree.get<int>("FiniteElementSetup.Parallel_Refine");

}

void Config::ReadAndInitMatProps()
{
    // Read MaterialProperties
    density = property_tree.get<double>("MaterialProperties.Density");
    Cp = property_tree.get<double>("MaterialProperties.Specific_Heat_Cp");
    CONDUCTIVITY_MODEL model = Conductivity_Model_Map.at(property_tree.get<string>("MaterialProperties.Thermal_Conductivity_Model"));
    switch (model)
    {
        case CONDUCTIVITY_MODEL::UNIFORM:
            //double kappa = property_tree.get<double>("MaterialProperties.Kappa");
            cond_model = new UniformCond(property_tree.get<double>("MaterialProperties.k"));
            break;
        /* TODO
        case CONDUCTIVITY_MODEL::LINEARIZED:
            cond_model = new LinearizedCond(property_tree.get<double>("MaterialProperties.k"), property_tree.get<double>("MaterialProperties.alpha"));
            break;
        */
    }
}

void Config::ReadIC()
{
    // Read InitialCondition
    BINARY_CHOICE restart_choice = Binary_Choice_Map.at(property_tree.get<string>("InitialCondition.Use_Restart"));
    use_restart = restart_choice == BINARY_CHOICE::YES ? true : false;
    // TODO: update to automatically retrieve from restarts, but for now:
    t0 = 0;

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
}
void Config::ReadpreCICE()
{
    // If no preCICE section detected, then do nothing + set with_preCICE
    if (property_tree.find("preCICE") == property_tree.not_found())
    {
        with_preCICE = false;
    }
    else // else preCICE section exists, read in info
    {
        with_preCICE = true;
        preCICE_participant_name = property_tree.get<string>("preCICE.Participant_Name");
        preCICE_config_file = property_tree.get<string>("preCICE.Config_File");
    }
}
void Config::ReadAndInitBCs(double& dt, BoundaryCondition** in_bcs, ParGridFunction* in_T_gf, SolverInterface* interface)
{
    // Read BoundaryConditions
    bc_count = property_tree.get_child("BoundaryConditions").size();
    in_bcs = new BoundaryCondition*[bc_count];
    
    // Assume no preCICE, update later if yes
    with_preCICE = false;
    
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
        //double value =  stod(bc_info[1].c_str());

        // Set the BC
        BOUNDARY_CONDITION type = Boundary_Condition_Map.at(bc_info[0]);
        
        double value;
        switch (type)
        {
            case BOUNDARY_CONDITION::HEATFLUX:
                value = stod(bc_info[1].c_str());
                in_bcs[index] = new UniformHeatFluxBC(attr, value);
                break;
            case BOUNDARY_CONDITION::ISOTHERMAL:
                value = stod(bc_info[1].c_str());
                in_bcs[index] = new UniformIsothermalBC(attr, value);
                break;
            case BOUNDARY_CONDITION::PRECICE_HEATFLUX:
                value = stod(bc_info[1].c_str());
                in_bcs[index] = new preCICEHeatFluxBC(attr, interface, in_T_gf, cond_model, dt, use_restart, preCICE_mesh_name, value);
                break;
            case BOUNDARY_CONDITION::PRECICE_ISOTHERMAL:
                value = stod(bc_info[1].c_str());
                in_bcs[index] = new preCICEIsothermalBC(attr, interface, in_T_gf, cond_model, dt, use_restart, preCICE_mesh_name, value);
                break;
        }

        index++;
        
        
        
    }
}

void Config::ReadTimeInt()
{
    // Read TimeIntegration
    time_scheme = Time_Scheme_Map.at(property_tree.get<string>("TimeIntegration.Time_Scheme"));
    dt = property_tree.get<double>("TimeIntegration.Delta_Time");
    tf = property_tree.get<double>("TimeIntegration.Final_Time");
}

void Config::ReadLinSolSettings()
{
    // Read LinearSolverSettings
    solver = Solver_Map.at(property_tree.get<string>("LinearSolverSettings.Solver"));
    prec = Preconditioner_Map.at(property_tree.get<string>("LinearSolverSettings.Preconditioner"));
    abs_tol = property_tree.get<double>("LinearSolverSettings.Absolute_Tolerance");
    rel_tol = property_tree.get<double>("LinearSolverSettings.Relative_Tolerance");
    max_iter = property_tree.get<int>("LinearSolverSettings.Max_Iterations");
}

void Config::ReadOutput()
{
    //Read Output
    restart_freq = property_tree.get<int>("Output.Restart_Freq");
    vis_freq = property_tree.get<int>("Output.Visualization_Freq");

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
            return "FMGRES";
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

Config::~Config()
{
    delete cond_model;
}