#pragma once
#include "material_property.hpp"

using namespace mfem;

class NonlinearJOTSDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
        ParFiniteElementSpace& fespace;
        MaterialProperty& k;
        ParGridFunction u_gf;
        GradientGridFunctionCoefficient grad_u_coeff;
        ScalarVectorProductCoefficient dkdu_times_grad_u;

        MixedScalarWeakDivergenceIntegrator term1;
        DiffusionIntegrator term2;


        Array<int> dofs;

    protected:
    public:
        NonlinearJOTSDiffusionIntegrator(ParFiniteElementSpace* fespace_, MaterialProperty& k_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
};


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
