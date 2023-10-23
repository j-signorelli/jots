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
        mfem::Coefficient* coeff;
    public:
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const

        virtual bool IsConstant() const = 0; // true if dk_du = 0
        virtual void UpdateCoeff(const mfem::Vector& T_ref) = 0;

        virtual std::string GetInitString() const = 0;

        virtual double GetLocalValue(double temp) const = 0;

        ~MaterialProperty() { delete coeff; };
};

class UniformProperty : public MaterialProperty
{
    private:
        const double k;
    protected:
    public:
        UniformProperty(const double in_k) : MaterialProperty(), k(in_k) { coeff = new mfem::ConstantCoefficient(k); };
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