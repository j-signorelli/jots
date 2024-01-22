#pragma once
#include "material_property.hpp"

using namespace mfem;

class JOTSNonlinearDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
    protected:
        Coefficient& lambda;
        Coefficient& dlambdadu;
        ParGridFunction u_gf;
        GradientGridFunctionCoefficient grad_u_coeff;
        ScalarVectorProductCoefficient dlambdadu_times_grad_u;

        MixedScalarWeakDivergenceIntegrator term1;
        DiffusionIntegrator term2;


        Array<int> dofs;
    public:
        JOTSNonlinearDiffusionIntegrator(ParFiniteElementSpace* fespace_, Coefficient& lambda_, Coefficient& dlambdadu_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
};

class JOTSNonlinearNeumannIntegrator : public NonlinearFormIntegrator
{
    private:
    protected:
        Coefficient& lambda;
        Coefficient& dlambdadu;

        BoundaryLFIntegrator vec_integ;
        BoundaryMassIntegrator grad_integ;
    
    public:
        JOTSNonlinearNeumannIntegrator(Coefficient& lambda_, Coefficient& dlambdadu_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
}

class JOTSNonlinearConvectionIntegrator : public NonlinearFormIntegrator
{
    private:
    protected:
        Coefficient& lambda;
        Coefficient& dlambdadu;
        ParGridFunction u_gf;
        GradientGridFunctionCoefficient grad_u_coeff;
        ScalarVectorProductCoefficient lambda_times_grad_u;
        InnerProductCoefficient dlambdadu_times_grad_u_dot_grad_u;

        MassIntegrator term1;
        ConvectionIntegrator term2;
    public:
        JOTSNonlinearConvectionIntegrator(ParFiniteElementSpace* fespace_, Coefficient& lambda_, Coefficient& dlambdadu_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
}