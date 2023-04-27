#pragma once

#include "mfem/mfem.hpp"


class SolverState
{
    protected:
        int it_num;
        double time;
        double dt;
        const double tf;
        mfem::Vector T;
        mfem::ParFiniteElementSpace& fespace;
        mutable mfem::ParGridFunction T_gf; // ONLY update whenever UpdateGF is called, instantiated once
        // mutable because we want to allow change to it in const fxns - it doesn't affect anyhting
    public:
        SolverState(mfem::ParGridFunction initial_gf, int in_it, double in_time, double in_dt, const double in_tf);

        int GetItNum() const { return it_num; };

        double GetTime() const { return time; };

        double& GetTimeRef() { return time; };

        double GetFinalTime() const { return tf; }

        double& GetdtRef() { return dt; };

        double Getdt() const { return dt; };

        mfem::ParFiniteElementSpace& GetParFESpace() const { return fespace; };

        mfem::Vector& GetTRef() { return T; };

        void UpdateGF() const { T_gf->SetFromTrueDofs(T); };

        mfem::ParGridFunction* GetGF() const { UpdateGF(); return T_gf; };

};