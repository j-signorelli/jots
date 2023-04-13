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
        double value;

    public:
        BoundaryCondition(int attr, double in_value, BOUNDARY_CONDITION in_type) : bdr_attr(attr), value(in_value), bc_type(in_type) {};
        int GetBdrAttr() const { return bdr_attr; }
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        double GetValue() const { return value; };
        virtual bool IsEssential() const = 0;
        virtual bool IsConstant() const = 0; // true if d/dt is 0 for this coefficient - required to ensure SetTime not called accidentally on something that doesn't have it
        virtual mfem::Coefficient* GetCoefficient() const= 0;
        virtual std::string GetInitString() const = 0;
};

class UniformIsothermalBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformIsothermalBC(int attr, double const_value) : BoundaryCondition(attr, const_value, BOUNDARY_CONDITION::ISOTHERMAL){};
        bool IsEssential() const { return true; }
        bool IsConstant() const { return true; }
        std::string GetInitString() const;
        mfem::Coefficient* GetCoefficient() const;
};

class UniformHeatFluxBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformHeatFluxBC(int attr, double const_value) : BoundaryCondition(attr, const_value, BOUNDARY_CONDITION::HEATFLUX){};
        bool IsEssential() const { return false; }
        bool IsConstant() const { return true; }
        std::string GetInitString() const;
        mfem::Coefficient* GetCoefficient() const;
};