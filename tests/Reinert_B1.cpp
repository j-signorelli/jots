#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;

double Reinert_B1_Analytical(const Vector& x, double time);

int main(int argc, char *argv[])
{  
    return Reinert_Test(1, Reinert_B1_Analytical);
}