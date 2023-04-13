#include "boundary_condition.hpp"

using namespace std;
using namespace mfem;

Coefficient* UniformIsothermalBC::GetCoefficient() const
{
    return new ConstantCoefficient(value);
}

string UniformIsothermalBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Isothermal --- Value: " << value;
    return sstm.str();
}

Coefficient* UniformHeatFluxBC::GetCoefficient() const
{
    return new ConstantCoefficient(value);
}

string UniformHeatFluxBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Heat Flux --- Value: " << value;
    return sstm.str();
}