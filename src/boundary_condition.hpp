#pragma once

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
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        double GetValue() const { return value; };
        virtual bool IsEssential() const = 0;
        virtual mfem::Coefficient* GetCoefficient() const= 0;

};

class UniformIsothermalBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformIsothermalBC(int attr, double const_value) : BoundaryCondition(attr, const_value, BOUNDARY_CONDITION::ISOTHERMAL){};
        bool IsEssential() const { return true; }
        mfem::Coefficient* GetCoefficient() const;
};

class UniformHeatFluxBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformHeatFluxBC(int attr, double const_value) : BoundaryCondition(attr, const_value, BOUNDARY_CONDITION::HEATFLUX){};
        bool IsEssential() const { return false; }
        mfem::Coefficient* GetCoefficient() const;
};