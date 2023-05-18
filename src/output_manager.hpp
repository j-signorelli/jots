#pragma once
#include <string>

#include "mfem.hpp"

#include "conductivity_model.hpp"

class OutputManager
{
    private:

        mfem::ConduitDataCollection* conduit_dc;
        mfem::ParaViewDataCollection* paraview_dc;
        mfem::ParGridFunction* rho_gf;
        mfem::ParGridFunction* Cp_gf;
        mfem::ParGridFunction* rank_gf;
        mfem::ParGridFunction* T_gf;
        mfem::ParGridFunction* k_gf;

        const mfem::Vector& T_ref;
        const ConductivityModel* cond_model; // Not allocated here

        void UpdateGridFunctions();

    protected:
    public:
        static const std::string TEMPERATURE;
        static const std::string RES_PREFIXPATH;
        static const std::string RES_NAME;
        
        OutputManager(mfem::ParFiniteElementSpace* fespace, const int fe_order, const double in_rho, const double in_Cp, const double in_rank, const mfem::Vector& in_T_ref, const ConductivityModel* in_cond_model);
        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);
        ~OutputManager();
};