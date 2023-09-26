include(FindMPI)
message(STATUS "Testing included!")
message(STATUS "Testing with MPI using: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")

configure_file(test_helper.hpp.in test_helper.hpp)

function(add_parallel_reg_test test_name test_main_file)
    add_executable(${test_name} ${test_main_file})
    target_link_libraries(${test_name} PUBLIC jotsLib)
    add_test(NAME ${test_name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ./${test_name})
endfunction()

function(add_precice_test dummy_name dummy_main_file jots_file_loc jots_file)

    # Create preCICE Dummy
    add_executable(${dummy_name} ${dummy_main_file})
    target_link_libraries(${dummy_name} PUBLIC jotsLib)
    add_test(NAME Reinert_B1_preCICE COMMAND ./${dummy_name} &
                                             (cd ${jots_file_loc} &&
                                                 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} jots -i ${jots_file}))

endfunction()

# Reinert_B1 tests "Transient Heat Conduction with Dirichlet Boundary Conditions"
add_parallel_reg_test(Reinert_B1 Reinert_B1.cpp)

# Reinert_B2 tests "Transient Heat Conduction with Neumann Boundary Conditions"
add_parallel_reg_test(Reinert_B2 Reinert_B2.cpp)

# Reinert_B3 tests "Transient Heat Conduction with Variable Thermal Properties"
add_parallel_reg_test(Reinert_B3 Reinert_B3.cpp)

# Restarted_Sinusoidal_Isothermal_Heated_Bar tests restart files
add_parallel_reg_test(Restart_Files Restart_Files.cpp)

# Reinert_B1_preCICE tests that temperature is received correctly in 1D
add_precice_test(Reinert_B1_preCICE_Dummy Reinert_B1_preCICE_Dummy.cpp ${PROJECT_SOURCE_DIR}/examples/Reinert_B1_preCICE/jots/ Reinert_B1_preCICE.ini)
