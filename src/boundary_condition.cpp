#include "boundary_condition.hpp"

using namespace std;
using namespace mfem;

string UniformIsothermalBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Isothermal --- Value: " << uniform_value;
    return sstm.str();
}

string UniformHeatFluxBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Heat Flux --- Value: " << uniform_value;
    return sstm.str();
}
/*
void UnsteadyNodalBC::InitCoefficient()
{
    // Get true boundary dofs
    fespace.GetBoundaryTrueDofs(boundary_dofs, attr);

    // Get boundary dof s
}
*/