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
        mfem::Coefficient* dcoeffdu;

    public:
        MaterialProperty() : coeff(nullptr), dcoeffdu(nullptr) {};
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const
        mfem::Coefficient& GetDCoeffRef() const { return *dcoeffdu; };
        virtual bool IsConstant() const = 0; // true if dcoeffdu = 0
        virtual void UpdateCoeff(const mfem::Vector& u_ref) = 0;
        virtual void UpdateCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs) = 0; // Both vectors same size; if none, then assume global
        virtual void UpdateDCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs) = 0;
        
        virtual std::string GetInitString() const = 0;
        virtual double GetLocalValue(double u_local) const = 0;
        virtual ~MaterialProperty() { delete coeff; delete dcoeffdu; };
};

class UniformProperty : public MaterialProperty
{
    private:
        const double mp_val;
    protected:
    public:
        UniformProperty(const double& in_mp);
        bool IsConstant() const { return true; }
        void UpdateCoeff(const mfem::Vector& u_ref) {};
        void UpdateCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs) {};
        void UpdateDCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs) {};
        std::string GetInitString() const;
        
        double GetLocalValue(double u_local) const { return mp_val; };
};

class PolynomialProperty : public MaterialProperty
{   
    private:
    protected:
        const std::vector<double> poly_coeffs;
        mfem::ParGridFunction mp_gf; // PGF for coeff
        mfem::ParGridFunction dmpdu_gf; // PGF for dcoeffdu
        
        mutable mfem::ParGridFunction u_gf;
        mutable mfem::ParGridFunction z;
    public:
        PolynomialProperty(const std::vector<double>& in_poly_coeffs, mfem::ParFiniteElementSpace& f);
        bool IsConstant() const { return false; };
        void UpdateCoeff(const mfem::Vector& u_ref);
        void UpdateCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs);
        void UpdateDCoeff(const mfem::Vector& u_ref_e, const mfem::Array<int>& dofs);

        std::string GetInitString() const;
        double GetLocalValue(double u_local) const;
};