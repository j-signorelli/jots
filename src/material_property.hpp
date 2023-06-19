#pragma once
#include <sstream>
#include <vector>
#include <cmath>

#include "mfem/mfem.hpp"

#include "option_structure.hpp"

class MaterialProperty
{   
    private:
    protected:
        MATERIAL_MODEL model;
        mfem::Coefficient* coeff;
    public:
        MaterialProperty(MATERIAL_MODEL in_model) : model(in_model) {}
        MATERIAL_MODEL GetModel() const { return model; };
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const

        virtual bool IsConstant() const = 0; // true if dk_dt = 0
        virtual void UpdateCoeff(const mfem::Vector& T_ref) = 0;

        virtual std::string GetInitString() const = 0;

        virtual double GetLocalValue(double temp) const = 0;

        ~MaterialProperty() { delete coeff; };
};

class UniformProperty : public MaterialProperty
{
    private:
        double k;
    protected:
    public:
        UniformProperty(double in_k) : MaterialProperty(MATERIAL_MODEL::UNIFORM), k(in_k) { coeff = new mfem::ConstantCoefficient(k); };
        bool IsConstant() const { return true; }
        void UpdateCoeff(const mfem::Vector& T_ref) {};

        std::string GetInitString() const;
        
        double GetLocalValue(double temp) const { return k; };
};

class PolynomialProperty : public MaterialProperty
{   
    private:
    protected:
        const std::vector<double> poly_coeffs;
        mfem::ParGridFunction* k_gf;

        mutable mfem::ParGridFunction T_gf;
        mutable mfem::ParGridFunction z;
    public:
        PolynomialProperty(const std::vector<double> in_poly_coeffs, mfem::ParFiniteElementSpace& f);
        bool IsConstant() const { return false; };
        void UpdateCoeff(const mfem::Vector& T_ref);

        std::string GetInitString() const;
        double GetLocalValue(double temp) const;
        ~PolynomialProperty() { delete k_gf; };
};