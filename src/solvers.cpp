#include "solvers.hpp"


bool UnsteadyHeatSolver::Running()
{
    return time < tf;
}