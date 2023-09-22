#include <cmath>

#include "jots_driver.hpp"

#include "test_helper.hpp"

using namespace std;
using namespace mfem;

double Reinert_B3_Analytical(const Vector& x, double time);

int main(int argc, char *argv[])
{  
    return Reinert_Test(3, Reinert_B3_Analytical);
}

double Reinert_B3_Analytical(const Vector& x, double time)
{
    double alpha = 2.5e-6;
    double L = 0.01;
    double k1 = 10.0;
    double k2 = 100;
    double T1 = 300;
    double T2 = 1300;

    double q_dot = 7.5e5;
    double T0 = 300;

    double D = ((k2-k1)/(T2-T1))*(1.0/(2.0*k1));
    double theta_0 = (T0-T1) + D*pow(T0-T1,2.0);

    double theta = alpha*time/pow(L,2.0) + 1.0/3.0 - x[0]/L + 0.5*pow(x[0]/L,2.0);

    for (int n = 1; n < N_REINERT+1; n++)
    {
        double A = (2.0/pow(M_PI,2.0))*(1.0/pow(n,2));
        double B = -pow(n,2.0)*pow(M_PI,2.0)*alpha*time/pow(L,2.0);
        double C = n*M_PI*x[0]/L;
        theta -= A*exp(B)*cos(C);
    }
    
    theta *= q_dot*L/k1;
    theta += theta_0;

    double a = D;
    double b = 1.0-2.0*D*T1;
    double c = D*pow(T1,2.0) - T1 - theta;

    return (-b+pow(pow(b,2.0)-4*a*c,0.5))/(2*a);

}