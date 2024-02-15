#pragma once
#include <string>
#include <map>

// This idea for arranging settings below was taken from SU2 v6.0:


enum class SIMULATION_TYPE
{
  LINEARIZED_UNSTEADY=0,
  NONLINEAR_UNSTEADY=1,
  STEADY=2 // Only Nonlinear Newton solver
};
static const std::map<std::string, SIMULATION_TYPE> Simulation_Type_Map = {{"Linearized_Unsteady", SIMULATION_TYPE::LINEARIZED_UNSTEADY},
                                                                           {"Nonlinear_Unsteady", SIMULATION_TYPE::NONLINEAR_UNSTEADY},
                                                                           {"Steady", SIMULATION_TYPE::STEADY}};


enum MATERIAL_PROPERTY : int
{
  DENSITY=0,
  SPECIFIC_HEAT=1,
  THERMAL_CONDUCTIVITY=2
};

static const std::map<std::string, MATERIAL_PROPERTY> Material_Property_Map = {{"Density_rho", MATERIAL_PROPERTY::DENSITY},
                                                                               {"Specific_Heat_C", MATERIAL_PROPERTY::SPECIFIC_HEAT},
                                                                               {"Thermal_Conductivity_k", MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY}};

enum class TIME_SCHEME
{

    // Implicit L-stable methods
    EULER_IMPLICIT = 0,
    //case 2:  ode_solver = new SDIRK23Solver(2); break;
    //case 3:  ode_solver = new SDIRK33Solver; break;

    // Explicit methods
    EULER_EXPLICIT = 1,
    RK4 = 2
    //case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
    //case 13: ode_solver = new RK3SSPSolver; break;
    //case 14: ode_solver = new RK4Solver; break;
    //case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;

    // Implicit A-stable methods (not L-stable)
    //case 22: ode_solver = new ImplicitMidpointSolver; break;
    //case 23: ode_solver = new SDIRK23Solver; break;
    //case 24: ode_solver = new SDIRK34Solver; break;
};

static const std::map<std::string, TIME_SCHEME> Time_Scheme_Map = {{"Euler_Implicit", TIME_SCHEME::EULER_IMPLICIT},
                                                         {"Euler_Explicit", TIME_SCHEME::EULER_EXPLICIT},
                                                         {"RK4", TIME_SCHEME::RK4}};

enum class BOUNDARY_CONDITION
{
  HEATFLUX,
  ISOTHERMAL,
  SINUSOIDAL_ISOTHERMAL,
  SINUSOIDAL_HEATFLUX,
  PRECICE_HEATFLUX,
  PRECICE_ISOTHERMAL
};

static const std::map<std::string, BOUNDARY_CONDITION> Boundary_Condition_Map = {{"HeatFlux", BOUNDARY_CONDITION::HEATFLUX},
                                                                       {"Isothermal", BOUNDARY_CONDITION::ISOTHERMAL},
                                                                       {"Sinusoidal_Isothermal", BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL},
                                                                       {"Sinusoidal_HeatFlux", BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX},
                                                                       {"preCICE_HeatFlux",  BOUNDARY_CONDITION::PRECICE_HEATFLUX},
                                                                       {"preCICE_Isothermal", BOUNDARY_CONDITION::PRECICE_ISOTHERMAL}};

enum class BINARY_CHOICE
{
  NO = 0,
  YES = 1
};

static const std::map<std::string, BINARY_CHOICE> Binary_Choice_Map = {{"Yes", BINARY_CHOICE::YES},
                                                             {"No", BINARY_CHOICE::NO}};

enum class MATERIAL_MODEL
{
  UNIFORM=0,
  POLYNOMIAL=1
};

static const std::map<std::string, MATERIAL_MODEL> Material_Model_Map = {{"Uniform", MATERIAL_MODEL::UNIFORM},
                                                                             {"Polynomial", MATERIAL_MODEL::POLYNOMIAL}};
                                                                    
enum class SOLVER
{
  CG=0,
  GMRES=1,
  FGMRES=2
};

static const std::map<std::string, SOLVER> Solver_Map = {{"CG", SOLVER::CG},
                                                         {"GMRES", SOLVER::GMRES},
                                                         {"FGMRES", SOLVER::FGMRES}};


enum class PRECONDITIONER
{
  JACOBI=0,
  CHEBYSHEV=16
};

static const std::map<std::string, PRECONDITIONER> Preconditioner_Map = {{"Jacobi", PRECONDITIONER::JACOBI},
                                                                         {"Chebyshev", PRECONDITIONER::CHEBYSHEV}};
