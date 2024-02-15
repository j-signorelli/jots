#include "../test_helper.hpp"

int main(int argc, char *argv[])
{  
    return Analytical_Reg_Test("NL_Reinert_B1.ini", Reinert_B1_Analytical, 1e-8); //h^p = (0.0001010101^2 ~ 1e-8)
}