#pragma once

#include "mfem.hpp"

#include "config_file.hpp"
#include "conductivity_model.hpp"

class OutputManager
{
    private:
        const static int RESTART_PRECISION;
        mfem::ParaViewDataCollection* paraview_dc;
        mfem::ParGridFunction* rho_gf;
        mfem::ParGridFunction* Cp_gf;
        mfem::ParGridFunction* rank_gf;
        mfem::ParGridFunction* T_gf;
        mfem::ParGridFunction* k_gf;

        const int rank;

        const mfem::Vector& T_ref;
        const ConductivityModel* cond_model; // Not allocated here
        const std::string output_restart_name;

        void UpdateGridFunctions();

    protected:
    public:
        OutputManager(const int in_rank, mfem::ParFiniteElementSpace* fespace, const Config* user_input, const mfem::Vector& in_T_ref, const ConductivityModel* in_cond_model);
        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);
        static std::tuple<double, double> GetTimeCyclesFromRestart(const std::string restart_info_line); // This must be consistent with how they are outputted

        ~OutputManager();
};