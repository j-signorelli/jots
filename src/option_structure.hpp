#pragma once
#include <string>
#include <map>

// This idea for arranging settings below was taken from SU2 v6.0.

// Because std::map::size() is not an integral-constant expression (evaluated at runtime), 
// although the sizes and values of each map below are known at compile-time,
// an alternative approach using the std::map iteration constructor is used for enum int types,
// where a size is needed for initialization.

// Pairs in each map are specified in an array
// See https://stackoverflow.com/questions/65670462/elegant-way-to-ensure-a-stdmap-has-a-concrete-size-in-compilation-time
// and https://stackoverflow.com/questions/2172053/c-can-i-statically-initialize-a-stdmap-at-compile-time


enum PHYSICS_TYPE : int
{
    THERMAL=0,
    STRUCTURAL=1
};

static const std::map<std::string, PHYSICS_TYPE>::value_type Physics_Type_Pairs[] = {{"Thermal", PHYSICS_TYPE::THERMAL},
                                                                                       {"Structural", PHYSICS_TYPE::STRUCTURAL}};
static const int PHYSICS_TYPE_SIZE = end(Physics_Type_Pairs) - begin(Physics_Type_Pairs);
static const std::map<std::string, PHYSICS_TYPE> Physics_Type_Map(begin(Physics_Type_Pairs), end(Physics_Type_Pairs));

static const std::map<PHYSICS_TYPE, std::string> Solution_Names_Map = {{PHYSICS_TYPE::THERMAL, "Temperature"},
                                                                        {PHYSICS_TYPE::STRUCTURAL, "Displacement"}};

enum class SIMULATION_TYPE
{
  LINEARIZED_UNSTEADY_HEAT,
  NONLINEAR_UNSTEADY_HEAT,
  STEADY_HEAT,
  EQUILIBRIUM_LINEAR_ELASTIC
};
static const std::map<std::string, SIMULATION_TYPE> Simulation_Type_Map = {{"Linearized_Unsteady_Heat", SIMULATION_TYPE::LINEARIZED_UNSTEADY_HEAT},
                                                                        {"Nonlinear_Unsteady_Heat", SIMULATION_TYPE::NONLINEAR_UNSTEADY_HEAT},
                                                                        {"Steady_Heat", SIMULATION_TYPE::STEADY_HEAT},
                                                                        {"Equilibrium_Linear_Elastic", SIMULATION_TYPE::EQUILIBRIUM_LINEAR_ELASTIC}};


enum MATERIAL_PROPERTY : int
{
  DENSITY=0,
  SPECIFIC_HEAT=1,
  THERMAL_CONDUCTIVITY=2,
  LAMES_FIRST_PARAMETER=3,
  LAMES_SECOND_PARAMETER=4
};
static const std::map<std::string, MATERIAL_PROPERTY>::value_type Material_Property_Pairs[] = {{"Density_rho", MATERIAL_PROPERTY::DENSITY},
                                                                                               {"Specific_Heat_C", MATERIAL_PROPERTY::SPECIFIC_HEAT},
                                                                                               {"Thermal_Conductivity_k", MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY},
                                                                                               {"Lames_First_Parameter_lambda", MATERIAL_PROPERTY::LAMES_FIRST_PARAMETER},
                                                                                               {"Lames_Second_Parameter_mu", MATERIAL_PROPERTY::LAMES_SECOND_PARAMETER}};
static const int MATERIAL_PROPERTY_SIZE = end(Material_Property_Pairs) - begin(Material_Property_Pairs);
static const std::map<std::string, MATERIAL_PROPERTY> Material_Property_Map(begin(Material_Property_Pairs), end(Material_Property_Pairs));



enum class TIME_SCHEME
{
    EULER_IMPLICIT,
    EULER_EXPLICIT,
    RK4
};

static const std::map<std::string, TIME_SCHEME> Time_Scheme_Map = {{"Euler_Implicit", TIME_SCHEME::EULER_IMPLICIT},
                                                                                 {"Euler_Explicit", TIME_SCHEME::EULER_EXPLICIT},
                                                                                 {"RK4", TIME_SCHEME::RK4}};

enum class THERMAL_BOUNDARY_CONDITION
{
  HEATFLUX,
  ISOTHERMAL,
  SINUSOIDAL_ISOTHERMAL,
  SINUSOIDAL_HEATFLUX,
  PRECICE_HEATFLUX,
  PRECICE_ISOTHERMAL
};

static const std::map<std::string, THERMAL_BOUNDARY_CONDITION> Thermal_Boundary_Condition_Map = {{"HeatFlux", THERMAL_BOUNDARY_CONDITION::HEATFLUX},
                                                                                            {"Isothermal", THERMAL_BOUNDARY_CONDITION::ISOTHERMAL},
                                                                                            {"Sinusoidal_Isothermal", THERMAL_BOUNDARY_CONDITION::SINUSOIDAL_ISOTHERMAL},
                                                                                            {"Sinusoidal_HeatFlux", THERMAL_BOUNDARY_CONDITION::SINUSOIDAL_HEATFLUX},
                                                                                            {"preCICE_HeatFlux",  THERMAL_BOUNDARY_CONDITION::PRECICE_HEATFLUX},
                                                                                            {"preCICE_Isothermal", THERMAL_BOUNDARY_CONDITION::PRECICE_ISOTHERMAL}};

enum class STRUCTURAL_BOUNDARY_CONDITION
{
    DISPLACEMENT,
    TRACTION,
    COMPONENT_WISE
};

static const std::map<std::string, STRUCTURAL_BOUNDARY_CONDITION> Structural_Boundary_Condition_Map = {{"Displacement", STRUCTURAL_BOUNDARY_CONDITION::DISPLACEMENT},
                                                                                                       {"Traction", STRUCTURAL_BOUNDARY_CONDITION::TRACTION},
                                                                                                       {"Component_Wise", STRUCTURAL_BOUNDARY_CONDITION::COMPONENT_WISE}};

enum class BINARY_CHOICE
{
  NO,
  YES
};

static const std::map<std::string, BINARY_CHOICE> Binary_Choice_Map = {{"Yes", BINARY_CHOICE::YES},
                                                             {"No", BINARY_CHOICE::NO}};

enum class MATERIAL_MODEL
{
  UNIFORM,
  POLYNOMIAL
};

static const std::map<std::string, MATERIAL_MODEL> Material_Model_Map = {{"Uniform", MATERIAL_MODEL::UNIFORM},
                                                                             {"Polynomial", MATERIAL_MODEL::POLYNOMIAL}};
                                                                    

enum class ADDITIONAL_SETTING
{
    PLANE_MODEL
};
static const std::map<ADDITIONAL_SETTING, std::string> Additional_Setting_String_Map = {{ADDITIONAL_SETTING::PLANE_MODEL, "2D_Plane_Model"}};

enum class ADDITIONAL_OPTION
{
    OFF,
    PLANE_STRESS,
    PLANE_STRAIN
};
static const std::map<ADDITIONAL_OPTION, std::string> Additional_Option_String_Map = {{ADDITIONAL_OPTION::OFF, "Off"},
                                                                                      {ADDITIONAL_OPTION::PLANE_STRESS, "Plane_Stress"},
                                                                                      {ADDITIONAL_OPTION::PLANE_STRAIN, "Plane_Strain"}};

enum class LINEAR_SOLVER
{
  CG,
  GMRES,
  FGMRES
};

static const std::map<std::string, LINEAR_SOLVER> Linear_Solver_Map = {{"CG", LINEAR_SOLVER::CG},
                                                         {"GMRES", LINEAR_SOLVER::GMRES},
                                                         {"FGMRES", LINEAR_SOLVER::FGMRES}};


enum class PRECONDITIONER
{
  JACOBI,
  CHEBYSHEV
};

static const std::map<std::string, PRECONDITIONER> Preconditioner_Map = {{"Jacobi", PRECONDITIONER::JACOBI},
                                                                         {"Chebyshev", PRECONDITIONER::CHEBYSHEV}};

enum class PRINT_LEVEL
{
  NONE,
  WARNINGS,
  ERRORS,
  ITERATIONS,
  FIRSTANDLAST,
  SUMMARY,
  ALL
};

static const std::map<std::string, PRINT_LEVEL> Print_Level_Map = {{"None", PRINT_LEVEL::NONE},
                                                                    {"Warnings", PRINT_LEVEL::WARNINGS},
                                                                    {"Errors", PRINT_LEVEL::ERRORS},
                                                                    {"Iterations", PRINT_LEVEL::ITERATIONS},
                                                                    {"FirstAndLast", PRINT_LEVEL::FIRSTANDLAST},
                                                                    {"Summary", PRINT_LEVEL::SUMMARY},
                                                                    {"All", PRINT_LEVEL::ALL}};
