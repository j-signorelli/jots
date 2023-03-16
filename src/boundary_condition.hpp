#pragma once

#include "option_structure.hpp"

// NOTE: This will soon become an abstract class with Dirichlet, Neumann, preCICE, etc.. children later on
class BoundaryCondition
{
    private:
        double value;
        BOUNDARY_CONDITION bc_type;

    public:
        BoundaryCondition(BOUNDARY_CONDITION type, double val);

        double GetValue() { return value; };
        BOUNDARY_CONDITION GetType() { return bc_type; }

    protected:
};