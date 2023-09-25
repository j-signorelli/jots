// This is a C++ version of Reinert_B3_preCICE dummy

#include "test_helper.hpp"

int main(int argc, char **argv)
{
    return Reinert_preCICE_Dummy(3, "Temperature", "Heat-Flux", 0.0, 7.5e5, Reinert_B3_Analytical);
}