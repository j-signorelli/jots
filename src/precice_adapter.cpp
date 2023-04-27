#include "precice_adapter.hpp"

using namespace std;
using namespace precice;

PreciceAdapter::PreciceAdapter(string in_part_name, string in_config, const int r, const int s)
: saved_state(nullptr),
  precice_bcs(nullptr),
  participant_name(in_part_name),
  config_file(in_config),
  rank(r),
  size(s)
{
    interface = new SolverInterface(participant_name, config_file, rank, size);
    dim = interface->getDimensions()
}

void PreciceAdapter::AddPreciceBCs(BoundaryCondition** in_bcs, vector<int> precice_bc_indices)
{
    num_bcs = precice_bc_indices.size()
    precice_bcs = new BoundaryCondition*[num_bcs];

    // Add BCs to list
    for (int i = 0; i < num_bcs; i++)
    {
        precice_bcs[i] = (*PreciceBC)in_bcs[precice_bc_indices[i]];

        // Set mesh ID
        int mesh_idid = interface->getMeshID(precice_bcs[i]->GetMeshName());
        precice_bcs[i]->SetMeshID(id);

        // Set mesh vertices
        interface->setMeshVertices(mesh_id, precice_bcs[i]->GetCoords());

        // Set read/write data IDs
        string read_data_name = precice_bcs[i]->GetReadDataName();
        string write_data_name = precice_bcs[i]->GetWriteDataName();
        precice_bcs[i]->SetReadDataID(interface->getDataID(read_data_name, mesh_id));
        precice_bcs[i]->SetWriteDataID(interface->getDataID(write_data_name, mesh_id));
    }

}



void PreciceAdapter::SaveOldState(const SolverState* state)
{
    // create new solver state object if it doesn't exist
    if (saved_state == nullptr)
    {
        saved_state = new SolverState(state->GetParFESpace());
        // Set the final time - this won't change throughout the simulation
        saved_state->SetFinalTime(state->SetFinalTime());
    }

    // Save all values of relevance
    saved_state->SetItNum(state->GetItNum());
    saved_state->SetTime(state->GetTime());
    saved_state->GetTRef() = state->GetTRef(); // This copies data

}
void ReloadOldState(SolverState* state) const
{

}

PreciceAdapter::~PreciceAdapter()
{
    delete saved_state;
    delete[] preCICE_BCs;
}