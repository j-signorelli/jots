#pragma once

#include "mfem/mfem.hpp"
#include "boundary_condition.hpp"
#include "option_structure.hpp"

class JOTSIterator
{
    private:
    protected:
        mfem::ParFiniteElementSpace& fespace;
        mfem::ParLinearForm b;// LinearForm representing Neumann BCs
        mfem::Vector b_vec; // Neumann BC Vector

        mfem::Array<int> ess_tdof_list; // list of essential true dofs

    public:
        JOTSIterator(mfem::ParFiniteElementSpace& f_, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, int bc_count);
        virtual void UpdateNeumann();

        virtual void Iterate(mfem::Vector& u) = 0;
        virtual void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) = 0;
        virtual ~JOTSIterator() = default;
};