#include "conductivity_model.hpp"

using namespace std;
using namespace mfem;

string UniformCond::GetInitString() const
{
    stringstream sstm;
    sstm << "Uniform --- k = " << k;
    return sstm.str();
}

string PolynomialCond::GetInitString() const
{
    stringstream sstm;
    sstm << "Polynomial --- k = ";
    int exp = 0;
    for (int i = 0; i < poly_coeffs.size(); i++)
    {   
        sstm << poly_coeffs[i];
        exp = poly_coeffs.size() - i - 1;
        if (exp == 1)
            sstm << "T + ";
        else if (exp != 0)
            sstm << "T^" << exp << " + ";
    }
    return sstm.str();
}

PolynomialCond::PolynomialCond(const std::vector<double> in_poly_coeffs, ParFiniteElementSpace& f) 
: ConductivityModel(CONDUCTIVITY_MODEL::POLYNOMIAL), 
  poly_coeffs(in_poly_coeffs)
{
    k_gf = new ParGridFunction(&f);
    
    coeff = new GridFunctionCoefficient(k_gf);
}

void PolynomialCond::UpdateCoeff(const mfem::Vector& T_ref)
{   
    // auxiliary PGFs
    ParGridFunction T_gf(k_gf->ParFESpace());
    ParGridFunction z(k_gf->ParFESpace());

    T_gf.SetFromTrueDofs(T_ref);
    
    *k_gf = 0;
    for (int i = 0; i < poly_coeffs.size(); i++)
    {
        z = 1;
        for (int j = 0; j < poly_coeffs.size() - i - 1; j++)
            z *= T_gf;
        k_gf->Add(poly_coeffs[i],z); // Update the GF associated with the coefficient
    }
}

double PolynomialCond::GetLocalConductivity(double temp) const
{
    double k = 0;
    for (int i = 0; i < poly_coeffs.size(); i++)
        k += poly_coeffs[i]*pow(temp, poly_coeffs.size() - i - 1);

    return k;
}