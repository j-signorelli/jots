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

        mfem::ParFiniteElementSpace& fespace;
        std::map<std::string, std::pair<mfem::Coefficient&, mfem::ParGridFunction*>> coeff_output_map;
        std::map<std::string, std::pair<const mfem::Vector&, mfem::ParGridFunction*>> vector_output_map;

        const int rank;

        void UpdateGridFunctions();

    protected:
    public:
        OutputManager(const int in_rank, mfem::ParFiniteElementSpace& f, const Config& user_input);
        
        void RegisterCoefficient(const std::string output_name, mfem::Coefficient& coeff);
        void RegisterSolutionVector(const std::string output_name, const mfem::Vector& vec);

        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);

        const mfem::ParGridFunction* GetT_gf();

        ~OutputManager();
};