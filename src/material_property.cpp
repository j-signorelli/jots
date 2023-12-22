#include "material_property.hpp"

using namespace std;
using namespace mfem;

UniformProperty::UniformProperty(const double& in_mp)
: MaterialProperty(),
  mp_val(in_mp)
{
    coeff = new ConstantCoefficient(mp_val);
    dcoeffdu = new ConstantCoefficient(0);
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
  u_gf(&f),
  z(&f)
{
    coeff = new GridFunctionCoefficient(&mp_gf);
    dcoeffdu = new GridFunctionCoefficient(&dmpdu_gf);
}

void PolynomialProperty::UpdateCoeff(const mfem::Vector& u_ref)
{   
    u_gf.SetFromTrueDofs(u_ref);
    
    mp_gf = 0;

    for (size_t i = 0; i < poly_coeffs.size(); i++)
    {
        z = 1;
        for (size_t j = 0; j < poly_coeffs.size() - i - 1; j++)
            z *= u_gf;
        mp_gf.Add(poly_coeffs[i],z); // Update the GF associated with the coefficient
    }
}

void PolynomialProperty::UpdateCoeff(const mfem::Vector& u_ref_e, const Array<int>& dofs)
{   
    mp_gf.SetSubVector(dofs, 0.0);

    Vector z_e(u_ref_e.Size());

    for (size_t i = 0; i < poly_coeffs.size(); i++)
    {
        z_e = 1;
        for (size_t j = 0; j < poly_coeffs.size() - i - 1; j++)
            z_e *= u_ref_e;
        mp_gf.AddElementVector(dofs, poly_coeffs[i], z_e); // Update the GF associated with the coefficient
    }
}

void PolynomialProperty::UpdateDCoeff(const mfem::Vector& u_ref_e, const Array<int>& dofs)
{   
    dmpdu_gf.SetSubVector(dofs, 0.0);

    Vector z_e(u_ref_e.Size());

    for (size_t i = 0; i < poly_coeffs.size() - 1; i++)
    {
        z_e = poly_coeffs.size() - 1 - i;
        for (size_t j = 0; j < poly_coeffs.size() - i - 2; j++)
            z_e *= u_ref_e;
        dmpdu_gf.AddElementVector(dofs, poly_coeffs[i], z_e); // Update the GF associated with the coefficient
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
