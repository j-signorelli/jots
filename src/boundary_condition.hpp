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

    public:
        BoundaryCondition(int attr, BOUNDARY_CONDITION in_type) : bdr_attr(attr), bc_type(in_type) {};
        int GetBdrAttr() const { return bdr_attr; }
        BOUNDARY_CONDITION GetType() const { return bc_type; };

        virtual bool IsEssential() const = 0;
        virtual bool IsConstant() const = 0; // true if d/dt is 0 for this coefficient - required to ensure SetTime not called accidentally on something that doesn't have it
        virtual mfem::Coefficient* GetCoefficient() const= 0;
        virtual std::string GetInitString() const = 0;
};

class UniformConstantBC : public BoundaryCondition
{   
    private:
    protected:
        double uniform_value;
    public:
        UniformConstantBC(int attr, BOUNDARY_CONDITION in_type, double in_value) : BoundaryCondition(attr, in_type), uniform_value(in_value) {};
        bool IsConstant() const { return true; }
        double GetValue() const { return uniform_value; }
        mfem::Coefficient* GetCoefficient() const { return new mfem::ConstantCoefficient(uniform_value); };
        
        virtual std::string GetInitString() const = 0;
};

class UniformIsothermalBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformIsothermalBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::ISOTHERMAL, const_value){};
        bool IsEssential() const { return true; }
        std::string GetInitString() const;
};

class UniformHeatFluxBC : public UniformConstantBC
{
    private:
    protected:
    public:
        UniformHeatFluxBC(int attr, double const_value) : UniformConstantBC(attr, BOUNDARY_CONDITION::HEATFLUX, const_value){};
        bool IsEssential() const { return false; }
        std::string GetInitString() const;
};

/*
class UnsteadyBC : public BoundaryCondition
{
    private:
    protected:
        virtual double TDF(mfem::Vector &x, double time) const = 0;
    
    public:
        UnsteadyBC(int attr, BOUNDARY_CONDITION in_type) : BoundaryCondition(attr, in_type) {};
        bool IsConstant() const { return false; };
        mfem::Coefficient* GetCoefficient() const {return new FunctionCoefficient(TDF); };
        
        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
}

class UnsteadyDiscreteBC : UnsteadyBC // Piecewise-Nodal Unsteady
{
    protected:
        std::
    
        double TDF(mfem::Vector &x, double time) const {};
    
    public:
        UnsteadyDiscreteBC(int attr, BOUNDARY_CONDITION in_type) : BoundaryCondition(attr, in_type) {};

        virtual bool IsEssential() const = 0;
        virtual std::string GetInitString() const = 0;
}
*/