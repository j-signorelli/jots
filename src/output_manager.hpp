#pragma once

#include "mfem.hpp"

#include "config_file.hpp"
#include "material_property.hpp"

class OutputManager
{
    private:
        const static int RESTART_PRECISION;
        mfem::VisItDataCollection* visit_dc;
        mfem::ParaViewDataCollection* paraview_dc;
        mfem::ParGridFunction* rho_gf;
        mfem::ParGridFunction* C_gf;
        mfem::ParGridFunction* rank_gf;
        mfem::ParGridFunction* T_gf;
        mfem::ParGridFunction* k_gf;

        mfem::Coefficient& k_coeff;
        mfem::Coefficient& C_coeff;

        const int rank;

        const mfem::Vector& T_ref;
        void UpdateGridFunctions();

    protected:
    public:
        OutputManager(const int in_rank, mfem::ParFiniteElementSpace* fespace, const Config& user_input, const mfem::Vector& in_T_ref, const MaterialProperty* rho_prop, const MaterialProperty* C_prop, const MaterialProperty* k_prop);
        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);

        const mfem::ParGridFunction* GetT_gf();

        ~OutputManager();
};