#include "boundary_condition.hpp"

using namespace mfem;

Coefficient* UniformIsothermalBC::GetCoefficient() const
{
    return new ConstantCoefficient(value);
}

Coefficient* UniformHeatFluxBC::GetCoefficient() const
{
    return new ConstantCoefficient(value);
}