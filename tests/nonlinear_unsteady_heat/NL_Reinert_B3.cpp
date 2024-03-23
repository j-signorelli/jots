#include "../test_helper.hpp"

int main(int argc, char *argv[])
{  
    return Analytical_Thermal_Reg_Test("NL_Reinert_B3.ini", Reinert_B3_Analytical, 1e-8);
}