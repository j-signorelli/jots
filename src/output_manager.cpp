#include "output_manager.hpp"

using namespace std;
using namespace mfem;

OutputManager::OutputManager(const int in_rank, ParFiniteElementSpace& f, const Config& user_input)
: visit_dc(nullptr),
  paraview_dc(nullptr),
  fespace(f),
  rank(in_rank),
  rank_coeff(in_rank)
{   
    //------------------------------------------------
    // Set up VisIt outputting (restarts)
    visit_dc = new VisItDataCollection(user_input.GetRestartPrefix(), fespace.GetParMesh());
    visit_dc->SetLevelsOfDetail(user_input.GetFEOrder());
    visit_dc->SetFormat(DataCollection::PARALLEL_FORMAT);
    visit_dc->SetPrecision(16);
    #ifdef MFEM_USE_ZLIB
        visit_dc->SetCompression(true);
    #endif 
    //------------------------------------------------
    // Set up ParaView outputting
    paraview_dc = new ParaViewDataCollection("ParaView", fespace.GetParMesh());
    paraview_dc->UseRestartMode(user_input.UsesRestart());
    paraview_dc->SetLevelsOfDetail(user_input.GetFEOrder());
    paraview_dc->SetDataFormat(VTKFormat::BINARY);
    paraview_dc->SetHighOrderOutput(true);

    //------------------------------------------------
    // Instantiate rank output
    // Rank:
    RegisterCoefficient("Rank", rank_coeff);

}

void OutputManager::RegisterCoefficient(const string output_name, Coefficient& coeff)
{

    // Create new CoefficientOutput + PGF for it
    coeff_output_map[output_name] = new CoefficientOutput(coeff, &fespace);

    // Project coefficient onto it
    coeff_output_map[output_name]->pgf->ProjectCoefficient(coeff);

    // Register
    paraview_dc->RegisterField(output_name, coeff_output_map[output_name]->pgf);
}

void OutputManager::RegisterSolutionVector(const string output_name, const Vector& vec)
{
    // Create new grid function for it and add Vector reference and GF pair to map
    vector_output_map[output_name] = new VectorOutput(vec, &fespace);

    // Set PGF from vector
    vector_output_map[output_name]->pgf->SetFromTrueDofs(vec);

    // Register to both ParaView + Restarts
    paraview_dc->RegisterField(output_name, vector_output_map[output_name]->pgf);
    visit_dc->RegisterField(output_name, vector_output_map[output_name]->pgf);
}

void OutputManager::UpdateGridFunctions()
{
    // Update all coefficient-driven GFs
    for (auto& x : coeff_output_map) // x.second is the value of a given {key,value} pair
        x.second->pgf->ProjectCoefficient(x.second->coeff_ref);

    // Update all vector-driven GFs
    for (auto& x : vector_output_map)
        x.second->pgf->SetFromTrueDofs(x.second->vector_ref);

}

const ParGridFunction* OutputManager::GetVectorPGF(string vec_label)
{
    UpdateGridFunctions();
    return vector_output_map[vec_label]->pgf;
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
    visit_dc->SetCycle(it_num);
    visit_dc->SetTime(time);
    visit_dc->Save();
}

OutputManager::~OutputManager()
{
    delete visit_dc;
    delete paraview_dc;
    for (auto& x : coeff_output_map)
        delete x.second;
    for (auto& x : vector_output_map)
        delete x.second;   
}