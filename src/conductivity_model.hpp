#pragma once
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
        virtual mfem::Coefficient* ApplyModel(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const = 0;

};

class ConstantCond : public ConductivityModel
{
    private:
        double k;
    protected:
    public:
        ConstantCond(double in_k) : ConductivityModel(CONDUCTIVITY_MODEL::CONSTANT), k(in_k) {};
        double Getk() const { return k; };
        mfem::Coefficient* ApplyModel(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const;
};

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
        mfem::Coefficient* ApplyModel(mfem::ParFiniteElementSpace* fespace, const mfem::Vector &u) const;

};