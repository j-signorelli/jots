#pragma once
#include <sstream>

#include "mfem/mfem.hpp"

#include "option_structure.hpp"

class ConductivityModel
{   
    private:
    protected:
        CONDUCTIVITY_MODEL model;
        mfem::Coefficient* coeff;
    public:
        ConductivityModel(CONDUCTIVITY_MODEL in_model) : model(in_model) {}
        CONDUCTIVITY_MODEL GetModel() const { return model; };
        mfem::Coefficient* GetCoeffPtr() const { return coeff; };

        virtual bool IsConstant() const = 0; // true if dk_dt = 0
        virtual void InitCoefficient() = 0;
        virtual void UpdateCoeff() = 0;

        virtual std::string GetInitString() const = 0;
        virtual double GetLocalConductivity(const mfem::ElementTransformation& transf, const mfem::IntegrationPoint& ip) const = 0;
        ~ConductivityModel() { delete coeff; };
};

class UniformCond : public ConductivityModel
{
    private:
        double k;
    protected:
    public:
        UniformCond(double in_k) : ConductivityModel(CONDUCTIVITY_MODEL::UNIFORM), k(in_k) {};
        bool IsConstant() const { return true; }
        void InitCoefficient() { coeff = new mfem::ConstantCoefficient(k); };
        void UpdateCoeff() {};

        std::string GetInitString() const;
        double GetLocalConductivity(const mfem::ElementTransformation& transf, const mfem::IntegrationPoint& ip) const { return k; };
};
/*
class LinearizedCond : public ConductivityModel
{
    private:
        double k;
        double alpha;
        As similarly done for sinusoidal BCs, send reference to solution vector in constructor.
    protected:
    public:
        LinearizedCond(double in_k, double in_alpha) : ConductivityModel(CONDUCTIVITY_MODEL::LINEARIZED), k(in_k), alpha(in_alpha) {};
        double Getk() const { return k; };
        double Getalpha() const { return alpha; };
        mfem::Coefficient* GetCoefficient(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const;
        bool IsConstant() const { return false; } //dk_dt != 0
};
*/