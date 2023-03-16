#pragma once
#include "option_structure.hpp"

class ConductivityModel
{   
    private:
    protected:
        CONDUCTIVITY_MODEL model;
    public:
        ConductivityModel(CONDUCTIVITY_MODEL in_model) : model(in_model) {}
        CONDUCTIVITY_MODEL GetModel() const { return model; };
        virtual void ApplyModel() const {};

};

class ConstantCond : public ConductivityModel
{
    private:
        double k;
    protected:
    public:
        ConstantCond(double in_k) : ConductivityModel(CONDUCTIVITY_MODEL::CONSTANT), k(in_k) {};
        double Getk() const { return k; };
        void ApplyModel() const;
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
        void ApplyModel() const;

};