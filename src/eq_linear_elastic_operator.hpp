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
        class LambdaStarCoefficient : public Coefficient
        {
            protected:
                Coefficient &lambda, &mu;
            public:
                LambdaStarCoefficient(const MaterialProperty& lambda_, const MaterialProperty& mu_)
                : lambda(lambda_.GetCoeffRef()), mu(mu_.GetDCoeffRef()) {};
                double Eval(ElementTransformation &T, const IntegrationPoint &ip);
        };
        
        ParBilinearForm a;
        HypreParMatrix *A_mat;
        LambdaStarCoefficient* lambda_star; // for plane stress
    public:
        EquilibriumLinearElasticOperator(ParFiniteElementSpace& f_, const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty &lambda_prop, MaterialProperty &mu_prop);
        
        void Iterate(Vector& u);

        // No updates mid run -- steady just completes inner iterations in Iterate and closes
        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};

        ~EquilibriumLinearElasticOperator() { delete lambda_star; delete A_mat;};
};