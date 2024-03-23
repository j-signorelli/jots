#include "../test_helper.hpp"

int main(int argc, char *argv[])
{  
    return Analytical_Thermal_Reg_Test("NL_Reinert_B2.ini", Reinert_B2_Analytical, 1e-8);
}