set(steady_tests
    Danish_Model1
    Danish_Model1_Neumann
)

add_parallel_reg_tests("${steady_tests}" Steady_Heat)

# Danish_Model1 tests nonlinear steady heat transfer with Dirichlet BCs
# Danish_Model1 tests nonlinear steady heat transfer with Neumann BC