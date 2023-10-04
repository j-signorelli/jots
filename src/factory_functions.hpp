
#pragma once
#include "mfem/mfem.hpp"

#include "option_structure.hpp"

namespace Factory
{

mfem::IterativeSolver* GetSolver(std::string solver_label, MPI_Comm comm_);

mfem::HypreSmoother::Type GetPrec(std::string prec_label);

}

#include "factory_functions.inl"