#pragma once

#include "mfem.hpp"

#include "config_file.hpp"
#include "material_property.hpp"
#include "jots_common.hpp"

struct CoefficientOutput
{
    mfem::Coefficient& coeff_ref; //ProjectCoefficient not const! Cannot have as const
    mfem::ParGridFunction pgf;

    CoefficientOutput(mfem::Coefficient &in_coeff, mfem::ParFiniteElementSpace &f) : coeff_ref(in_coeff), pgf(&f) {};
};


struct VectorOutput
{
    const mfem::Vector& vector_ref;
    mfem::ParGridFunction pgf;

    VectorOutput(const mfem::Vector &in_vec, mfem::ParFiniteElementSpace &f) : vector_ref(in_vec), pgf(&f) {};
};

class OutputManager
{
    private:
        mfem::VisItDataCollection visit_dc;
        mfem::ParaViewDataCollection paraview_dc;

        // For visualizing parallel decomposition:
        mfem::L2_FECollection pw_const_fec;
        mfem::ParFiniteElementSpace pw_const_fes;
        mfem::ParGridFunction rank_pgf;

        // Scalar FESpace for Coefficient visualization
        mfem::ParFiniteElementSpace& scalar_fes;

        std::map<std::string, CoefficientOutput*> coeff_output_map;
        std::map<std::string, VectorOutput*> vector_output_map;

        void UpdateGridFunctions();

    protected:
    public:
        OutputManager(const int in_rank, const Config &user_input, mfem::ParFiniteElementSpace& in_scalar_fes);
        
        void RegisterCoefficient(const std::string output_name, mfem::Coefficient& coeff);
        void RegisterSolutionVector(const std::string output_name, const mfem::Vector& vec, mfem::ParFiniteElementSpace& f);

        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);

        const mfem::ParGridFunction* GetVectorPGF(std::string vec_label);

        ~OutputManager();
};