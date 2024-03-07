#pragma once
#include "jots_iterator.hpp"
#include "jots_common.hpp"
#include "material_property.hpp"
#include "config_file.hpp"

using namespace mfem;

class EquilibriumLinearElasticOperator : public JOTSIterator
{
    private:
    protected:
        ParBilinearForm a;
        HypreParMatrix *A_mat;
    public:
        EquilibriumLinearElasticOperator(ParFiniteElementSpace& f_, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty &lambda_prop, MaterialProperty &mu_prop);
        
        void Iterate(Vector& u);

        // No updates mid run -- steady just completes inner iterations in Iterate and closes
        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};

        ~EquilibriumLinearElasticOperator() { delete A_mat; };
};