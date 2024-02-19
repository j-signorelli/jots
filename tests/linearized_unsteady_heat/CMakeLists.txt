set(linearized_unsteady_tests
    Reinert_B1
    Reinert_B2
    Reinert_B3
    Restart_Files
    Reinert_B1_preCICE
    Reinert_B3_preCICE
)

add_parallel_reg_tests("${linearized_unsteady_tests}" Linearized_Unsteady_Heat)

# Reinert_B1 tests "Transient Heat Conduction with Dirichlet Boundary Conditions"
# Reinert_B2 tests "Transient Heat Conduction with Neumann Boundary Conditions"
# Reinert_B3 tests "Transient Heat Conduction with Variable Thermal Properties"
# Restarted_Sinusoidal_Isothermal_Heated_Bar tests restart files
# Reinert_B1_preCICE tests that HF is received correctly in 1D + temperature sent
# Reinert_B3_preCICE tests that temperature is received correctly in 1D + HF sent