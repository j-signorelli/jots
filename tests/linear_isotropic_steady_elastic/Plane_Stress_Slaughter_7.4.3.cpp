#include "../test_helper.hpp"

void Analytical_Slaughter_7p4p3(const Vector &x, double time, Vector &u)
{   
    double lambda = 1;
    double mu = 1;
    double L = 1.0;
    double h = 0.1;
    double P_per_w = 1e-3*h;
    double E = mu*(3*lambda+2*mu)/(lambda+mu);
    double nu = lambda/(2*(lambda+mu));
    double I_3_per_w = (1.0/12.0)*pow(h,3);

    double X_1 = x[0];
    double X_2 = x[1];

    u[0] = (-P_per_w/(2*E*I_3_per_w))*(pow(L,2) - pow(X_1,2))*X_2
           - ((P_per_w*(2+nu))/(6*E*I_3_per_w))*pow(X_2,3)
           + ((P_per_w*(1+nu)*pow(h,2))/(8*E*I_3_per_w))*X_2;
    u[1] = ((-P_per_w*pow(L,3))/(6*E*I_3_per_w))*(
           2 - (3*X_1/L)*(1-nu*pow(X_2,2)/pow(L,2))
           + pow(X_1,3)/pow(L,3) + ((3*pow(h,2))/(4*pow(L,2)))*(1+nu)*(1-X_1/L));
}

int main(int argc, char *argv[])
{  
    return Analytical_Structural_Reg_Test("Plane_Stress_Slaughter_7.4.3.ini", Analytical_Slaughter_7p4p3, 1e-4); //h^p = (0.01^2 = 1e-4)
}