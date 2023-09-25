// This is a C++ version of Reinert_B1_preCICE dummy

#include "test_helper.hpp"

int main(int argc, char **argv)
{
    return Reinert_preCICE_Dummy(1, "Heat-Flux", "Temperature", 0.01, 500, Reinert_B1_Analytical_HF);
}