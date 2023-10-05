#include "output_manager.hpp"

using namespace std;
using namespace mfem;

//const int OutputManager::RESTART_PRECISION = 16;

OutputManager::OutputManager(const int in_rank, ParFiniteElementSpace& f, const Config& user_input)
: rank(in_rank),
  fespace(f)
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
    ConstantCoefficient rank_coeff(in_rank);
    RegisterCoefficient("Rank", rank_coeff);

}

void OutputManager::RegisterCoefficient(const string output_name, Coefficient& coeff)
{

    // Create new grid function for it and add Coefficient reference and GF pair to map
    pair<Coefficient&, ParGridFunction*> coeff_pair(coeff, new ParGridFunction(&fespace));
    coeff_output_map[output_name] = coeff_pair;

    // Project coefficient onto it
    coeff_output_map[output_name].second->ProjectCoefficient(coeff);

    // Register
    paraview_dc->RegisterField(output_name, coeff_output_map[output_name].second);
}

void OutputManager::RegisterSolutionVector(const string output_name, const Vector& vec)
{
    // Create new grid function for it and add Vector reference and GF pair to map
    pair<const Vector&, ParGridFunction*> vec_pair(vec, new ParGridFunction(&fespace));
    vector_output_map[output_name] = vec_pair;

    // Project coefficient onto it
    vector_output_map[output_name].second->SetFromTrueDofs(vec);

    // Register to both ParaView + Restarts
    paraview_dc->RegisterField(output_name, vector_output_map[output_name].second);
    visit_dc->RegisterField(output_name, vector_output_map[output_name].second);
}

void OutputManager::UpdateGridFunctions()
{
    // Update all coefficient-driven GFs
    for (map<string, pair<Coefficient&, ParGridFunction*>> it = coeff_output_map.begin(); it != coeff_output_map; it++)
    {
        string key = it->first;
        coeff_output_map[key].second->ProjectCoefficient(coeff_output_map[key].first);
    }

    // Update all vector-driven GFs
    for (map<string, pair<Vector&, ParGridFunction*>> it = vector_output_map.begin(); it != vector_output_map; it++)
    {
        string key = it->first;
        vector_output_map[key].second->SetFromTrueDofs(vector_output_map[key].first);
    }
    
}

const ParGridFunction* OutputManager::GetT_gf()
{
    UpdateGridFunctions();
    return T_gf;
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
    for (map<string, pair<Coefficient&, ParGridFunction*>> it = coeff_output_map.begin(); it != coeff_output_map.end(); it++)
        delete it->second.second;
    for (map<string, pair<Vector&, ParGridFunction*>> it = vector_output_map.begin(); it != vector_output_map.end(); it++)
        delete it->second.second;   
}