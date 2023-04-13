#pragma once
#include <sstream>

#include "mfem/mfem.hpp"

#include "option_structure.hpp"

class ConductivityModel
{   
    private:
    protected:
        CONDUCTIVITY_MODEL model;
    public:
        ConductivityModel(CONDUCTIVITY_MODEL in_model) : model(in_model) {}
        CONDUCTIVITY_MODEL GetModel() const { return model; };
        virtual bool IsConstant() const = 0; // true if dk_dt = 0
        virtual std::string GetInitString() const = 0;
        virtual mfem::Coefficient* GetCoefficient() const = 0;//(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const = 0;
        
};

class UniformCond : public ConductivityModel
{
    private:
        double k;
    protected:
    public:
        UniformCond(double in_k) : ConductivityModel(CONDUCTIVITY_MODEL::UNIFORM), k(in_k) {};
        double Getk() const { return k; };
        bool IsConstant() const { return true; }
        std::string GetInitString() const;
        mfem::Coefficient*  GetCoefficient() const;
};
/* TODO:
class LinearizedCond : public ConductivityModel
{
    private:
        double k;
        double alpha;
    protected:
    public:
        LinearizedCond(double in_k, double in_alpha) : ConductivityModel(CONDUCTIVITY_MODEL::LINEARIZED), k(in_k), alpha(in_alpha) {};
        double Getk() const { return k; };
        double Getalpha() const { return alpha; };
        mfem::Coefficient* GetCoefficient(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const;
        bool IsConstant() const { return false; } //dk_dt != 0
};
*/