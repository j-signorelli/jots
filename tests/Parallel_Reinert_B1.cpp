#include <cmath>

#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;

double Reinert_B1_Analytical(const Vector& x, double time);

int main(int argc, char *argv[])
{  
    return Reinert_Test(1, Reinert_B1_Analytical);
}

double Reinert_B1_Analytical(const Vector& x, double time)
{
    double theta = 1.0;
    double alpha = 2.5e-6;
    double L = 0.01;

    for (int n = 0; n < N_REINERT+1; n++)
    {
        double A = pow(-1.0,n)/(2.0*n+1.0);
        double B = (-pow(2.0*n+1.0,2.0)*pow(M_PI,2)*alpha*time)/(4.0*pow(L,2));
        double C = ((2.0*n+1.0)*M_PI*x[0])/(2.0*L);

        theta = theta - (4.0/M_PI)*A*exp(B)*cos(C);
    }
    
    double T_0 = 300.0;
    double T_D = 500.0;

    return theta*(T_D - T_0) + T_0;

}