#pragma once
#include "jots_iterator.hpp"
#include "material_property.hpp"
using namespace mfem;

class NonlinearJOTSDiffusionIntegrator : public NonlinearFormIntegrator
{
    private:
        const MaterialProperty& k;
        ParFiniteElementSpace& fespace;

        ParGridFunction f;
        ParGridFunction df;
        GridFunctionCoefficient f_coeff;
        GridFunctionCoefficient df_coeff;

        DomainLFGradIntegrator integ_vec;
        DiffusionIntegrator ??? integ_vec;

        Array<int> dofs;
    protected:
    public:
        NonlinearJOTSDiffusionIntegrator(const MaterialProperty& k_, ParFiniteElementSpace* fespace_);
        void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect);
        void AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat);
}

class SteadyConductionOperator : public JOTSIterator
{
    private:
    protected:
        const MaterialProperty& k;
        NewtonSolver
    public:
        SteadyConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, const MaterialProperty& k_prop);
        void Iterate(mfem::Vector& u);


        // No updates mid run -- steady just completes inner iterations in Iterate and closes
        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) = {};
        void UpdateNeumann() = {};
}