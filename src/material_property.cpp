#include "material_property.hpp"

using namespace std;
using namespace mfem;

string UniformProperty::GetInitString() const
{
    stringstream sstm;
    sstm << "Uniform --- Value = " << k;
    return sstm.str();
}

string PolynomialProperty::GetInitString() const
{
    stringstream sstm;
    sstm << "Polynomial --- Value = ";
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

PolynomialProperty::PolynomialProperty(const std::vector<double> in_poly_coeffs, ParFiniteElementSpace& f) 
: MaterialProperty(), 
  poly_coeffs(in_poly_coeffs),
  T_gf(&f),
  z(&f)
{
    k_gf = new ParGridFunction(&f);
    
    coeff = new GridFunctionCoefficient(k_gf);
}

void PolynomialProperty::UpdateCoeff(const mfem::Vector& T_ref)
{   
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

double PolynomialProperty::GetLocalValue(double temp) const
{
    double k = 0;
    for (int i = 0; i < poly_coeffs.size(); i++)
        k += poly_coeffs[i]*pow(temp, poly_coeffs.size() - i - 1);

    return k;
}