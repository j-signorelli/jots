
#pragma once
#include "mfem/mfem.hpp"

#include "option_structure.hpp"

namespace Factory
{

mfem::IterativeSolver* GetSolver(std::string solver_label, MPI_Comm comm_);

mfem::HypreSmoother::Type GetPrec(std::string prec_label);

mfem::ODESolver* GetODESolver(std::string time_scheme_label);

}

template<typename Key, typename Value>
std::vector<Key> GetKeyVector(std::map<Key, Value> in_map);

#include "helper_functions.inl"