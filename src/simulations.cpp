#include "simulations.hpp"



bool UnsteadyHeatSimulation::UsesMaterialProperty(MATERIAL_PROPERTY mat_prop)
{
    switch (mat_prop)
    {
        case MATERIAL_PROPERTY::DENSITY:
            return true;
            break;
        case MATERIAL_PROPERTY::SPECIFIC_HEAT:
            return true;
            break;
        case MATERIAL_PROPERTY::THERMAL_CONDUCTIVITY:
            return true;
            break;
        default:
            return false;
    }
}