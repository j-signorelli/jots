#pragma once

#include "mfem.hpp"

#include "option_structure.hpp"

class BoundaryCondition
{
    private:
    protected:
        BOUNDARY_CONDITION bc_type;
        double value; // TEMPORARY

    public:
        BoundaryCondition(double in_value, BOUNDARY_CONDITION in_type) : value(in_value), bc_type(in_type) {};
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        double GetValue() const { return value; }; // TEMPORARY
        virtual bool IsEssential() const = 0;
        virtual mfem::Coefficient* GetCoefficient() const= 0;

};

class UniformIsothermalBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformIsothermalBC(double const_value) : BoundaryCondition(const_value, BOUNDARY_CONDITION::ISOTHERMAL){};
        bool IsEssential() const { return true; }
        mfem::Coefficient* GetCoefficient() const;
};

class UniformHeatFluxBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformHeatFluxBC(double const_value) : BoundaryCondition(const_value, BOUNDARY_CONDITION::HEATFLUX){};
        bool IsEssential() const { return false; }
        mfem::Coefficient* GetCoefficient() const;
};