#include "../test_helper.hpp"

int main(int argc, char *argv[])
{  
    return Analytical_Thermal_Reg_Test("Danish_Model1.ini", Danish_Model1_Analytical,3e-5);
}

