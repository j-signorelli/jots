#pragma once
#include "material_property.hpp"

using namespace mfem;

class JOTSNonlinearDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
        ParFiniteElementSpace& fespace;
        Coefficient& lambda;
        Coefficient& dlambdadu;
        ParGridFunction u_gf;
        GradientGridFunctionCoefficient grad_u_coeff;
        ScalarVectorProductCoefficient dlambdadu_times_grad_u;

        MixedScalarWeakDivergenceIntegrator term1;
        DiffusionIntegrator term2;


        Array<int> dofs;

    protected:
    public:
        JOTSNonlinearDiffusionIntegrator(ParFiniteElementSpace* fespace_, Coefficient& lambda_, Coefficient& dlambdadu_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
};

/*
class NonlinearJOTSMassIntegrator : public NonlinearFormIntegrator
{
    private:
        ParFiniteElementSpace& fespace;
        MaterialProperty& rho;
        MaterialProperty& C;
        ParGridFunction u_gf;
        GridFunctionCoefficient u_coeff;
        ProductCoefficient drhodu_C;
        ProductCoefficient rho_dCdu;
        SumCoefficient mat_prop_coeff;
        ProductCoefficient mat_prop_coeff_times_u;

        MassIntegrator term1;
        MassIntegrator term2;

        Array<int> dofs;

    protected:
    public:
        NonlinearJOTSMassIntegrator(ParFiniteElementSpace* fespace_, MaterialProperty& rho_, MaterialProperty& C_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
    
};
*/