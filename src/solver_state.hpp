#pragma once

#include "mfem/mfem.hpp"


class SolverState
{
    protected:
        int it_num;
        double time;
        double dt;
        double tf;
        mfem::Vector T;
        mfem::ParFiniteElementSpace* fespace; // Not allocated here
        mutable mfem::ParGridFunction* T_gf; // ONLY update whenever UpdateGF is called, instantiated once
        // mutable because we want to allow change to it in const fxns - it doesn't affect anyhting
    public:
        SolverState(mfem::ParFiniteElementSpace* f);

        void SetItNum(int in_it) { it_num = in_it; };

        int GetItNum() const { return it_num; };

        void SetTime(double in_time) { time = in_time; };

        double GetTime() const { return time; };

        double& GetTimeRef() { return time; };

        void SetFinalTime(double in_tf) { tf = in_tf; };

        double GetFinalTime() const { return tf; }

        void Setdt(double in_dt) { dt = in_dt; };

        double& GetdtRef() { return dt; };

        double Getdt() const { return dt; };

        mfem::ParFiniteElementSpace* GetParFESpace() const { return fespace; }; // Note that we are allowing FESpace to be changed externally here despite this fxn being declared constant, but won't for JOTS

        mfem::Vector& GetTRef() { return T; };

        void UpdateGF() const { T_gf->SetFromTrueDofs(T); };

        mfem::ParGridFunction* GetGF() const { UpdateGF(); return T_gf; };

        ~SolverState();

};