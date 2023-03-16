#include "boundary_condition.hpp"

BoundaryCondition::BoundaryCondition(BOUNDARY_CONDITION type, double val)
 :  value(val),
    bc_type(type)

{
    // Nothing else (as of now) needed
    // May want to pass additional stuff later, but might be OK to just implement functions that take pointers to shit
}