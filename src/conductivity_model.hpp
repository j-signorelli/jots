#pragma once
#include <sstream>
#include <vector>
#include <cmath>

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
        virtual void UpdateCoeff(const mfem::Vector& T_ref) = 0;

        virtual std::string GetInitString() const = 0;

        virtual double GetLocalConductivity(double temp) const = 0;

        ~ConductivityModel() { delete coeff; };
};

class UniformCond : public ConductivityModel
{
    private:
        double k;
    protected:
    public:
        UniformCond(double in_k) : ConductivityModel(CONDUCTIVITY_MODEL::UNIFORM), k(in_k) { coeff = new mfem::ConstantCoefficient(k); };
        bool IsConstant() const { return true; }
        void UpdateCoeff(const mfem::Vector& T_ref) {};

        std::string GetInitString() const;
        
        double GetLocalConductivity(double temp) const { return k; };
};

class PolynomialCond : public ConductivityModel
{   
    private:
    protected:
        const std::vector<double> poly_coeffs;
        mfem::ParGridFunction* T_gf;
        mfem::ParGridFunction* k_gf;

    public:
        PolynomialCond(const std::vector<double> in_poly_coeffs, mfem::ParFiniteElementSpace& f);
        bool IsConstant() const { return false; };
        void UpdateCoeff(const mfem::Vector& T_ref);

        std::string GetInitString() const;
        double GetLocalConductivity(double temp) const;
        ~PolynomialCond() { delete T_gf; delete k_gf; };
};