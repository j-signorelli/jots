#include "conductivity_model.hpp"

using namespace mfem;

Coefficient* UniformCond::GetCoefficient() const
{
    return new ConstantCoefficient(k);
}

/*
Coefficient* LinearizedCond::GetCoefficient(ParFiniteElementSpace* fespace, const Vector &u) const
{
    ParGridFunction u_alpha_gf(fespace);
    u_alpha_gf.SetFromTrueDofs(u);
    for (int i = 0; i < u_alpha_gf.Size(); i++)
    {
        u_alpha_gf(i) = k + alpha*u_alpha_gf(i);
    }

    return new GridFunctionCoefficient(&u_alpha_gf);
}
*/