#pragma once
#include "material_property.hpp"

using namespace mfem;

class NonlinearJOTSDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
        MaterialProperty& k;
        ParFiniteElementSpace& fespace;
        ParGridFunction u_gf;
        GradientGridFunctionCoefficient grad_u_coeff;
        ScalarVectorProductCoefficient dkdu_times_grad_u;

        MixedScalarWeakDivergenceIntegrator term1;
        DiffusionIntegrator term2;


        Array<int> dofs;

    protected:
    public:
        NonlinearJOTSDiffusionIntegrator(MaterialProperty& k_, ParFiniteElementSpace* fespace_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
};

/*
class NonlinearJOTSMassIntegrator : public NonlinearFormIntegrator
{
    private:
        MaterialProperty& rho;
        MaterialProperty& C;
        ParFiniteElementSpace& fespace;
        ParGridFunction u_gf;
        GridFunctionCoefficient u_coeff;
        ProductCoefficient rho_C;
        ProductCoefficient rho_dCdu;
        ProductCoefficient rho_dCdu_times_u;

        MassIntegrator mass;
        MassIntegrator mass_d;

        Array<int> dofs;

    protected:
    public:
        NonlinearJOTSMassIntegrator(MaterialProperty& rho_, MaterialProperty& C_, ParFiniteElementSpace* fespace_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
    
};
*/