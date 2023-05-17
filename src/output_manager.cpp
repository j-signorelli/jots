#include "output_manager.hpp"

using namespace mfem;

OutputManager::OutputManager(mfem::ParFiniteElementSpace* fespace, const int fe_order, const double in_rho, const double in_Cp, const double in_rank, const Vector& in_T_ref, const ConductivityModel* in_cond_model)
: T_ref(in_T_ref),
  cond_model(in_cond_model)
{   
    //------------------------------------------------
    // Set up restart file outputting

    
    //------------------------------------------------
    // Set up ParaView outputting
    paraview_dc = new ParaViewDataCollection("Output", fespace->GetParMesh());
    // TODO: Flag as restart or not -- should be able to just append to an existing pvd file.
    paraview_dc->SetPrefixPath("ParaView");
    paraview_dc->SetLevelsOfDetail(fe_order);
    paraview_dc->SetDataFormat(VTKFormat::BINARY);
    paraview_dc->SetHighOrderOutput(true);

    //------------------------------------------------
    // Instantiate output grid functions + register:
    // Rank:
    rank_gf = new ParGridFunction(fespace);
    ConstantCoefficient rank_coeff(in_rank);
    rank_gf->ProjectCoefficient(rank_coeff);
    paraview_dc->RegisterField("Rank", rank_gf);

    // Density:
    rho_gf = new ParGridFunction(fespace);
    ConstantCoefficient rho_coeff(in_rho);
    rho_gf->ProjectCoefficient(rho_coeff);
    paraview_dc->RegisterField("Density", rho_gf);

    // Specific Heat:
    Cp_gf = new ParGridFunction(fespace);
    ConstantCoefficient Cp_coeff(in_Cp);
    Cp_gf->ProjectCoefficient(Cp_coeff);
    paraview_dc->RegisterField("Specific_Heat", Cp_gf);

    // Thermal Conductivity:
    k_gf = new ParGridFunction(fespace);
    k_gf->ProjectCoefficient(*cond_model->GetCoeffPtr());
    paraview_dc->RegisterField("Thermal_Conductivity", k_gf);

    // Temperature:
    T_gf = new ParGridFunction(fespace);
    T_gf->SetFromTrueDofs(T_ref);
    paraview_dc->RegisterField("Temperature", T_gf);

    //------------------------------------------------
    

}

void OutputManager::UpdateGridFunctions()
{
    // Update temperature GF right off the bat
    T_gf->SetFromTrueDofs(T_ref);

    // If non-constant thermal conductivity, update k GF with updated coefficient
    if (!cond_model->IsConstant())
        k_gf->ProjectCoefficient(*cond_model->GetCoeffPtr());

}

void OutputManager::WriteVizOutput(const int it_num, const double time)
{   
    UpdateGridFunctions();
    paraview_dc->SetCycle(it_num);
    paraview_dc->SetTime(time);
    paraview_dc->Save();
}

OutputManager::~OutputManager()
{
    delete paraview_dc;
    delete rho_gf;
    delete Cp_gf;
    delete rank_gf;
    delete T_gf;
    delete k_gf;
}