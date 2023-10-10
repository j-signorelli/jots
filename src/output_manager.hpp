#pragma once

#include "mfem.hpp"

#include "config_file.hpp"
#include "material_property.hpp"
#include "helper_functions.hpp"

struct CoefficientOutput
{
    mfem::Coefficient& coeff_ref; //ProjectCoefficient not const! Cannot have as const
    mfem::ParGridFunction* pgf;

    CoefficientOutput(mfem::Coefficient& in_coeff, mfem::ParFiniteElementSpace* f) : coeff_ref(in_coeff), pgf(nullptr) { pgf = new mfem::ParGridFunction(f); };
    ~CoefficientOutput() { delete pgf; };
};


struct VectorOutput
{
    const mfem::Vector& vector_ref;
    mfem::ParGridFunction* pgf;

    VectorOutput(const mfem::Vector& in_vec, mfem::ParFiniteElementSpace* f) : vector_ref(in_vec), pgf(nullptr) { pgf = new mfem::ParGridFunction(f); };
    ~VectorOutput() { delete pgf; };
};

class OutputManager
{
    private:
        mfem::VisItDataCollection* visit_dc;
        mfem::ParaViewDataCollection* paraview_dc;

        mfem::ParFiniteElementSpace& fespace;
        mfem::ConstantCoefficient rank_coeff;

        std::map<std::string, CoefficientOutput*> coeff_output_map;
        std::map<std::string, VectorOutput*> vector_output_map;

        const int rank;

        void UpdateGridFunctions();

    protected:
    public:
        OutputManager(const int in_rank, mfem::ParFiniteElementSpace& f, const Config& user_input);
        
        void RegisterCoefficient(const std::string output_name, mfem::Coefficient& coeff);
        void RegisterSolutionVector(const std::string output_name, const mfem::Vector& vec);

        void WriteVizOutput(const int it_num, const double time);
        void WriteRestartOutput(const int it_num, const double time);

        const mfem::ParGridFunction* GetVectorPGF(std::string vec_label);

        ~OutputManager();
};