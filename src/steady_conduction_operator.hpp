#pragma once
#include "jots_iterator.hpp"
#include "material_property.hpp"
#include "config_file.hpp"

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

class SteadyConductionOperator : public JOTSIterator
{
    private:
    protected:
        ParNonlinearForm k;

        IterativeSolver* lin_solver;
        HypreSmoother lin_prec;
        NewtonSolver newton;
        
    public:
        SteadyConductionOperator(const Config& in_config, const BoundaryCondition* const* in_bcs, mfem::Array<int>* all_bdr_attr_markers, MaterialProperty& k_prop, ParFiniteElementSpace& f_);
        void Iterate(mfem::Vector& u);


        // No updates mid run -- steady just completes inner iterations in Iterate and closes
        void ProcessMatPropUpdate(MATERIAL_PROPERTY mp) {};
        void UpdateNeumann() {};

        ~SteadyConductionOperator() { delete lin_solver; };
};