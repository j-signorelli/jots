#include "conductivity_model.hpp"

using namespace std;
using namespace mfem;

string UniformCond::GetInitString() const
{
    stringstream sstm;
    sstm << "Uniform --- k: " << k;
    return sstm.str();
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