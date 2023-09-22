#include <cmath>

#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;

double Reinert_B2_Analytical(const Vector& x, double time);

int main(int argc, char *argv[])
{  
    return Reinert_Test(2, Reinert_B2_Analytical);
}

double Reinert_B2_Analytical(const Vector& x, double time)
{
    double alpha = 2.5e-6;
    double L = 0.01;
    double k = 10.0;
    double q_dot = 7.5e5;

    double theta = alpha*time/pow(L,2.0) + 1.0/3.0 - x[0]/L + 0.5*pow(x[0]/L,2.0);

    for (int n = 1; n < N_REINERT+1; n++)
    {
        double A = -pow(n,2.0)*pow(M_PI,2.0)*alpha*time/pow(L,2.0);
        double B = n*M_PI*x[0]/L;
        double C = (2.0/pow(M_PI,2.0))*(1.0/pow(n,2.0));
        theta -= C*exp(A)*cos(B);
    }
    
    double T_0 = 300.0;

    return theta*q_dot*L/k + T_0;

}