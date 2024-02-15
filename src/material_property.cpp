#include "material_property.hpp"

using namespace std;
using namespace mfem;

UniformProperty::UniformProperty(const double& in_mp)
: MaterialProperty(),
  mp_val(in_mp)
{
    coeff = new ConstantCoefficient(mp_val);
    dcoeffdu = new ConstantCoefficient(0);
    d2coeffdu2 = new ConstantCoefficient(0);
}

string UniformProperty::GetInitString() const
{
    stringstream sstm;
    sstm << "Uniform --- Value = " << mp_val;
    return sstm.str();
}

PolynomialProperty::PolynomialProperty(const std::vector<double>& in_poly_coeffs, ParFiniteElementSpace& f) 
: MaterialProperty(), 
  poly_coeffs(in_poly_coeffs),
  mp_gf(&f),
  dmpdu_gf(&f),
  d2mpdu2_gf(&f),
  z1(f.GetTrueVSize()),
  z2(f.GetTrueVSize())
{
    coeff = new GridFunctionCoefficient(&mp_gf);
    dcoeffdu = new GridFunctionCoefficient(&dmpdu_gf);
    d2coeffdu2 = new GridFunctionCoefficient(&d2mpdu2_gf);
}

void PolynomialProperty::UpdateCoeff(const mfem::Vector& u_ref)
{   
    z1 = 0.0;
    for (size_t i = 0; i < poly_coeffs.size(); i++)
    {
        z2 = 1.0;
        for (size_t j = 0; j < poly_coeffs.size() - i - 1; j++)
            z2 *= u_ref;
        z2 *= poly_coeffs[i];
        z1 += z2;
    }

    mp_gf.SetFromTrueDofs(z1); // tdofs --> dofs
}

void PolynomialProperty::UpdateDCoeff(const mfem::Vector& u_ref)
{   
    z1 = 0.0;
    for (size_t i = 0; i < poly_coeffs.size() - 1; i++)
    {
        z2 = poly_coeffs.size() - 1 - i;
        for (size_t j = 0; j < poly_coeffs.size() - i - 2; j++)
            z2 *= u_ref;
        z2 *= poly_coeffs[i];
        z1 += z2;
    }

    dmpdu_gf.SetFromTrueDofs(z1);
}

void PolynomialProperty::UpdateD2Coeff(const mfem::Vector& u_ref)
{   
    z1 = 0.0;

    for (size_t i = 0; i < poly_coeffs.size() - 2; i++)
    {
        z2 = (poly_coeffs.size() - 1 - i) * (poly_coeffs.size() - 2 - i);
        for (size_t j = 0; j < poly_coeffs.size() - i - 3; j++)
            z2 *= u_ref;
        z2 *= poly_coeffs[i];
        z1 += z2;
        
    }

    d2mpdu2_gf.SetFromTrueDofs(z1);
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
