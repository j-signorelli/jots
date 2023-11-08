#include "material_property.hpp"

using namespace std;
using namespace mfem;

UniformProperty::UniformProperty(const double& in_mp)
: MaterialProperty(),
  mp_val(in_mp)
{
    coeff = new ConstantCoefficient(mp_val);
}

string UniformProperty::GetInitString() const
{
    stringstream sstm;
    sstm << "Uniform --- Value = " << k;
    return sstm.str();
}

PolynomialProperty::PolynomialProperty(const std::vector<double>& in_poly_coeffs, ParFiniteElementSpace& f) 
: MaterialProperty(), 
  poly_coeffs(in_poly_coeffs),
  mp_gf(&f),
  u_gf(&f),
  z(&f)
{
    coeff = new GridFunctionCoefficient(mp_gf);
}

void PolynomialProperty::UpdateCoeff(const mfem::Vector& u_ref)
{   
    u_gf.SetFromTrueDofs(u_ref);
    
    *mp_gf = 0;
    for (size_t i = 0; i < poly_coeffs.size(); i++)
    {
        z = 1;
        for (size_t j = 0; j < poly_coeffs.size() - i - 1; j++)
            z *= u_gf;
        mp_gf->Add(poly_coeffs[i],z); // Update the GF associated with the coefficient
    }
}

string PolynomialProperty::GetInitString() const
{
    stringstream sstm;
    sstm << "Polynomial --- Value = ";
    int exp = 0;
    for (size_t i = 0; i < poly_coeffs.size(); i++)
    {   
        sstm << poly_coeffs[i];
        exp = poly_coeffs.size() - i - 1;
        if (exp == 1)
            sstm << "u + ";
        else if (exp != 0)
            sstm << "u^" << exp << " + ";
    }
    return sstm.str();
}


double PolynomialProperty::GetLocalValue(double u_local) const
{
    double k = 0;
    for (size_t i = 0; i < poly_coeffs.size(); i++)
        k += poly_coeffs[i]*pow(u_local, poly_coeffs.size() - i - 1);

    return k;
}