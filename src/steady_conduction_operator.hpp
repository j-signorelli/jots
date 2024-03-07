#pragma once
#include "jots_iterator.hpp"
#include "jots_common.hpp"
#include "material_property.hpp"
#include "config_file.hpp"
#include "jots_nlfis.hpp"

using namespace mfem;

class SteadyConductionOperator : public JOTSIterator
{
    private:
    protected:
        ParNonlinearForm k;
        JOTSNewtonSolver newton;
        
    public:
        SteadyConductionOperator(ParFiniteElementSpace& f_, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& k_prop);
        void Iterate(Vector& u);


        // No updates mid run -- steady just completes inner iterations in Iterate and closes
        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};
};