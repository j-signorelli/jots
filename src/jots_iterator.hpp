#pragma once

#include "mfem/mfem.hpp"
#include "boundary_condition.hpp"
#include "option_structure.hpp"
#include "config_file.hpp"
#include "jots_common.hpp"

class JOTSIterator
{
    private:
    protected:
        mfem::ParFiniteElementSpace& fespace;
        mfem::ParLinearForm b;// LinearForm representing Neumann BC assuming user-input value = coefficient
        mfem::Vector b_vec; // Neumann BC Vector (tdofs)
        mfem::VectorArrayCoefficient neumann_coeff;
        mfem::Array<int> ess_tdof_list; // list of essential true dofs

        mfem::IterativeSolver *lin_solver;    // Linear solver
        mfem::HypreSmoother lin_prec; // Preconditioner for linear solver
        JOTSNewtonSolver newton;

    public:
        JOTSIterator(mfem::ParFiniteElementSpace& f_, const Config &in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, int bc_count);
        virtual void UpdateNeumann();

        virtual void Iterate(mfem::Vector& u) = 0;
        virtual void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) = 0;
        virtual ~JOTSIterator() { delete lin_solver; }; // Explicitly define destructor as virtual, so derived types can have their own destructors
};