include(FindMPI)
message(STATUS "Testing included!")
message(STATUS "Testing with MPI using: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}")

configure_file(test_helper.hpp.in test_helper.hpp)

function(add_parallel_reg_test test_name test_category)
    configure_file(${test_name}.ini.in ${test_name}.ini @ONLY)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} PUBLIC jotsLib)
    add_test(NAME ${test_category}/${test_name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ./${test_name})
    set_tests_properties(${test_category}/${test_name} PROPERTIES LABELS ${test_category})
endfunction()

function(add_parallel_reg_tests tests tests_category)
    foreach(test ${tests})
        add_parallel_reg_test(${test} ${tests_category})
    endforeach()
endfunction()

add_subdirectory(linearized_unsteady_heat)
add_subdirectory(nonlinear_unsteady_heat)
add_subdirectory(steady_heat)