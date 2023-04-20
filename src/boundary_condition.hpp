#pragma once
#include <sstream>

#include "mfem.hpp"

#include "option_structure.hpp"

class BoundaryCondition
{
    private:
    protected:
        int bdr_attr;
        BOUNDARY_CONDITION bc_type;
        mfem::Coefficient* coeff;
        
    public:
        BoundaryCondition(int attr, BOUNDARY_CONDITION in_type) : bdr_attr(attr), bc_type(in_type) {};
        int GetBdrAttr() const { return bdr_attr; }
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        mfem::Coefficient* GetCoeffPtr() const { return coeff; };

        virtual void InitCoefficient() = 0;
        virtual void UpdateCoeff() = 0;
        virtual bool IsEssential() const = 0;
        virtual bool IsConstant() const = 0; // true if d/dt is 0 for this coefficient
        virtual std::string GetInitString() const = 0;

        ~BoundaryCondition() {delete coeff;};
};

class UniformConstantBC : public BoundaryCondition
{   
    private:
    protected:
        double uniform_value;
    public:
        UniformConstantBC(int attr, BOUNDARY_CONDITION in_type, double in_value) : BoundaryCondition(attr, in_type), uniform_value(in_value) {};
        bool IsConstant() const { return true; };
        double GetValue() const { return uniform_value; };
        void InitCoefficient() { coeff = new mfem::ConstantCoefficient(uniform_value); };
        void UpdateCoeff() {};
        virtual std::string GetInitString() const = 0;
};

class UniformIsothermalBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformIsothermalBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::ISOTHERMAL, const_value){};
        bool IsEssential() const { return true; };
        std::string GetInitString() const;
};

class UniformHeatFluxBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformHeatFluxBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::HEATFLUX, const_value){};
        bool IsEssential() const { return false; };
        std::string GetInitString() const;
};

/*
class UnsteadyNodalBC : public BoundaryCondition
{
    private:
    protected:
        mfem::ParFiniteElementSpace &fespace;
        mfem::ParGridFunction coeff_gf;
        mfem::Array<int> boundary_dofs; // From FiniteElementSpace::GetBoundaryTrueDofs
        mfem::Array<double> 

    public:
        UnsteadyNodalBC(int attr, BOUNDARY_CONDITION in_type, mfem::ParFiniteElementSpace &f) : BoundaryCondition(attr, in_type), fespace(f), coeff_gf(f) {};
        bool IsConstant() const { return false; };
        void InitCoefficient();
        
        void UpdateCoeff();

        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
};
*/