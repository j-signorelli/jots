#include "output_manager.hpp"

using namespace std;
using namespace mfem;

const int OutputManager::RESTART_PRECISION = 15;

OutputManager::OutputManager(const int in_rank, ParFiniteElementSpace* fespace, const Config* user_input, const Vector& in_T_ref, const ConductivityModel* in_cond_model)
: rank(in_rank),
  T_ref(in_T_ref),
  cond_model(in_cond_model),
  output_restart_name(user_input->GetOutputRestartFile())
{   
    //------------------------------------------------
    // Set up ParaView outputting
    paraview_dc = new ParaViewDataCollection("ParaView", fespace->GetParMesh());
    paraview_dc->UseRestartMode(user_input->UsesRestart());
    paraview_dc->SetLevelsOfDetail(user_input->GetFEOrder());
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
    ConstantCoefficient rho_coeff(user_input->GetDensity());
    rho_gf->ProjectCoefficient(rho_coeff);
    paraview_dc->RegisterField("Density", rho_gf);

    // Specific Heat:
    Cp_gf = new ParGridFunction(fespace);
    ConstantCoefficient Cp_coeff(user_input->GetCp());
    Cp_gf->ProjectCoefficient(Cp_coeff);
    paraview_dc->RegisterField("Specific_Heat", Cp_gf);

    // Thermal Conductivity:
    k_gf = new ParGridFunction(fespace);
    k_gf->ProjectCoefficient(*cond_model->GetCoeffPtr());
    paraview_dc->RegisterField("Thermal_Conductivity", k_gf);

    // Temperature
    //------------------------------------------------
    T_gf = new ParGridFunction(fespace);
    T_gf->SetFromTrueDofs(T_ref);
    paraview_dc->RegisterField("Temperature", T_gf);

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

void OutputManager::WriteRestartOutput(const int it_num, const double time)
{
    UpdateGridFunctions();

    ofstream out_file;
    out_file << fixed << setprecision(RESTART_PRECISION);

    stringstream sstm;
    sstm << output_restart_name << "_" << it_num << ".dat";
    out_file.open(sstm.str(), ios::out);
    out_file << time << endl;
    out_file << it_num << endl;
    // Only save T_gf to restart
    T_gf->SaveAsOne(out_file);
    
    out_file.close();
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