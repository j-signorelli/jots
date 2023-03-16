#pragma once

#include "option_structure.hpp"

class BoundaryCondition
{
    private:
    protected:
        BOUNDARY_CONDITION bc_type;
        double value;

    public:
        BoundaryCondition(double in_value, BOUNDARY_CONDITION in_type) : value(in_value), bc_type(in_type) {};
        BOUNDARY_CONDITION GetType() const { return bc_type; };
        double GetValue() { return value; };
        virtual void ApplyBC() const {};

};


class UniformIsothermalBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformIsothermalBC(double const_value) : BoundaryCondition(const_value, BOUNDARY_CONDITION::ISOTHERMAL){};
        void ApplyBC() const;
};

class UniformHeatFluxBC : public BoundaryCondition
{
    private:
    protected:
    public:
        UniformHeatFluxBC(double const_value) : BoundaryCondition(const_value, BOUNDARY_CONDITION::HEATFLUX){};
        void ApplyBC() const;
};