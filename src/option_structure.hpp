#pragma once
#include <string>
#include <map>

// This idea for arranging settings below was taken from SU2 v6.0:

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
  PRECICE_HEATFLUX,
  PRECICE_ISOTHERMAL
};

static const std::map<std::string, BOUNDARY_CONDITION> Boundary_Condition_Map = {{"HeatFlux", BOUNDARY_CONDITION::HEATFLUX},
                                                                       {"Isothermal", BOUNDARY_CONDITION::ISOTHERMAL},
                                                                       {"Sinusoidal_Isothermal", BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL},
                                                                       {"preCICE_HeatFlux",  BOUNDARY_CONDITION::PRECICE_HEATFLUX},
                                                                       {"preCICE_Isothermal", BOUNDARY_CONDITION::PRECICE_ISOTHERMAL}};

enum class BINARY_CHOICE
{
  NO = 0,
  YES = 1
};

static const std::map<std::string, BINARY_CHOICE> Binary_Choice_Map = {{"Yes", BINARY_CHOICE::YES},
                                                             {"No", BINARY_CHOICE::NO}};

enum class CONDUCTIVITY_MODEL
{
  UNIFORM=0,
  LINEARIZED=1
};

static const std::map<std::string, CONDUCTIVITY_MODEL> Conductivity_Model_Map = {{"Uniform", CONDUCTIVITY_MODEL::UNIFORM},
                                                                             {"Linearized", CONDUCTIVITY_MODEL::LINEARIZED}};
                                                                    
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
