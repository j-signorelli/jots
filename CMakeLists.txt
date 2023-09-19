# Require same CMake as MFEM
cmake_minimum_required(VERSION 3.10)

# Include CTest for regression tests
#include(CTest)

# Define our project
project(jots VERSION 1.1 LANGUAGES CXX)

# Create interface library for compiler flag requirements (requiring C++11)
add_library(compiler_flag_requirements INTERFACE)
target_compile_features(compiler_flag_requirements INTERFACE cxx_std_11)

# Get preCICE
find_package(PRECICE REQUIRED CONFIG NAMES precice)
message(STATUS "Found preCICE config: TRUE")

# Get MFEM
find_package(MFEM REQUIRED CONFIG)
message(STATUS "Found MFEM config: TRUE")

# Get Boost headers
find_package(Boost REQUIRED CONFIG COMPONENTS headers)
message(STATUS "Found Boost headers: TRUE")

# Add src directory (creates jotsLib)
add_subdirectory(src)

# Create executable
add_executable(jots jots_main.cpp)
target_link_libraries(jots PUBLIC jotsLib)

# Install
install(TARGETS jots DESTINATION bin)
install(FILES "${PROJECT_SOURCE_DIR}/pw_mfem_mesh_fix.py" DESTINATION bin)

# Add regression tests:
# Serial_Reinert_B1 tests "Transient Heat Conduction with Dirichlet Boundary Conditions"
#add_test(NAME Serial_Reinert_B1 COMMAND )

# Reinert_B2 tests "Transient Heat Conduction with Neumann Boundary Conditions"
#add_test(NAME Reinert_B2 COMMAND )

# Reinert_B3 tests "Transient Heat Conduction with Variable Thermal Properties"
#add_test(NAME Reinert_B3 COMMAND )

# Modal_Oscillation_Plate_preCICE tests specification of isothermal BC from preCICE
#add_test(NAME Modal_Oscillation_Plate_preCICE COMMAND )

# Restarted_Sinusoidal_Isothermal_Heated_Bar tests restart files
#add_test(NAME Restarted_Sinusoidal_Isothermal_Heated_Bar COMMAND )