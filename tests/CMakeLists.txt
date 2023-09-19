include(FindMPI)
message(STATUS "Testing included!")
message(STATUS "Testing with MPI using: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4")

# Parallel_Reinert_B1 tests "Transient Heat Conduction with Dirichlet Boundary Conditions"
add_test(NAME Parallel_Reinert_B1 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 jots -i ${CMAKE_SOURCE_DIR}/examples/Reinert_B1/Reinert_B1.ini)

# Reinert_B2 tests "Transient Heat Conduction with Neumann Boundary Conditions"
#add_test(NAME Reinert_B2 COMMAND )

# Reinert_B3 tests "Transient Heat Conduction with Variable Thermal Properties"
#add_test(NAME Reinert_B3 COMMAND )

# Modal_Oscillation_Plate_preCICE tests specification of isothermal BC from preCICE
#add_test(NAME Modal_Oscillation_Plate_preCICE COMMAND )

# Restarted_Sinusoidal_Isothermal_Heated_Bar tests restart files
#add_test(NAME Restarted_Sinusoidal_Isothermal_Heated_Bar COMMAND )