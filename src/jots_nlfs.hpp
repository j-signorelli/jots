#pragma once
#include "material_property.hpp"

using namespace mfem;

class NonlinearJOTSDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
        MaterialProperty& k;
        ParFiniteElementSpace& fespace;
        ParGridFunction u_gf;
        GridFunctionCoefficient u_coeff;
        ProductCoefficient dkdu_times_u;

        DiffusionIntegrator diff;
        DiffusionIntegrator diff_d;


        Array<int> dofs;

    protected:
    public:
        NonlinearJOTSDiffusionIntegrator(MaterialProperty& k_, ParFiniteElementSpace* fespace_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
};