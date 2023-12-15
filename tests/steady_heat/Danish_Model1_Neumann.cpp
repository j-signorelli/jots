#include "../test_helper.hpp"

int main(int argc, char *argv[])
{  
    return Analytical_Reg_Test("Danish_Model1_Neumann.ini", Danish_Model1_Analytical, 3e-5);
}
