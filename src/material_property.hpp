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
        mfem::Coefficient* d2coeffdu2;
    public:
        MaterialProperty() : coeff(nullptr), dcoeffdu(nullptr), d2coeffdu2(nullptr) {};
        mfem::Coefficient& GetCoeffRef() const { return *coeff; }; // To be used only for assigning to linear/bilinear forms or projecting coeff, so declared const
        mfem::Coefficient& GetDCoeffRef() const { return *dcoeffdu; };
        mfem::Coefficient& GetD2CoeffRef() const { return *d2coeffdu2; };

        void UpdateAllCoeffs(const mfem::Vector& u_ref) { UpdateCoeff(u); UpdateDCoeff(u); UpdateD2Coeff(u); };

        virtual bool IsConstant() const = 0; // true if dcoeffdu = 0
        virtual void UpdateCoeff(const mfem::Vector& u_ref) = 0;
        virtual void UpdateDCoeff(const mfem::Vector& u_ref) = 0;
        virtual void UpdateD2Coeff(const mfem::Vector& u_ref) = 0;

        virtual std::string GetInitString() const = 0;
        virtual double GetLocalValue(double u_local) const = 0;
        virtual ~MaterialProperty() { delete coeff; delete dcoeffdu; delete d2coeffdu; };
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
        void UpdateDCoeff(const mfem::Vector& u_ref) {};
        void UpdateD2Coeff(const mfem::Vector& u_ref) {};
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
        mfem::ParGridFunction d2mpdu2_gf; // PGF for d2coeffdu2

        mutable mfem::ParGridFunction z;
    public:
        PolynomialProperty(const std::vector<double>& in_poly_coeffs, mfem::ParFiniteElementSpace& f);
        bool IsConstant() const { return false; };
        void UpdateCoeff(const mfem::Vector& u_ref);
        void UpdateDCoeff(const mfem::Vector& u_ref);
        void UpdateD2Coeff(const mfem::Vector& u_ref);

        std::string GetInitString() const;
        double GetLocalValue(double u_local) const;
};