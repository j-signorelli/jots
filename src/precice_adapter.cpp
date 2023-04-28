#include "precice_adapter.hpp"

using namespace std;
using namespace precice;
using namespace mfem;

const string PreciceAdapter::cowid = precice::constants::actionWriteInitialData();
const string PreciceAdapter::cowic = precice::constants::actionWriteIterationCheckpoint();
const string PreciceAdapter::coric = precice::constants::actionReadIterationCheckpoint();

PreciceAdapter::PreciceAdapter(string in_part_name, string in_config, const int r, const int s)
: precice_bcs(nullptr),
  participant_name(in_part_name),
  config_file(in_config),
  rank(r),
  size(s)
{
    interface = new SolverInterface(participant_name, config_file, rank, size);
    dim = interface->getDimensions();
}

void PreciceAdapter::AddPreciceBCs(BoundaryCondition** in_bcs, vector<int> precice_bc_indices)
{
    num_bcs = precice_bc_indices.size();
    precice_bcs = new PreciceBC*[num_bcs];

    // Add BCs to list
    for (int i = 0; i < num_bcs; i++)
    {
        precice_bcs[i] = dynamic_cast<PreciceBC*>(in_bcs[precice_bc_indices[i]]);

        // Set mesh ID
        int mesh_id = interface->getMeshID(precice_bcs[i]->mesh_name);
        precice_bcs[i]->mesh_id = mesh_id;

        // Set mesh vertices
        interface->setMeshVertices(mesh_id, precice_bcs[i]->num_dofs, precice_bcs[i]->coords, precice_bcs[i]->vertex_ids);

        // Set read/write data IDs
        precice_bcs[i]->read_data_id = interface->getDataID(precice_bcs[i]->read_data_name, mesh_id);
        precice_bcs[i]->write_data_id = interface->getDataID(precice_bcs[i]->write_data_name, mesh_id);
    }

}

void PreciceAdapter::WriteInitialData(const Vector T, const ConductivityModel* cond_model)
{
    // For isothermal wall, sending heat flux
    //      If restart, get heatflux from state, send as is.
    //      If not restart, must project coeff first, then get heat flux, then send

    // For heatflux wall, sending temperature
    //      If restart, get temperature from state and send
    //      If not restart, not at all any difference.
    for (int i = 0; i < num_bcs; i++)
    {
        precice_bcs[i]->RetrieveInitialWriteData(T, cond_model);

        interface->writeBlockScalarData(precice_bcs[i]->write_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->write_data_arr);
    }
}

void PreciceAdapter::GetReadData()
{
    for (int i = 0; i < num_bcs; i++)
    {
        interface->readBlockScalarData(precice_bcs[i]->read_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->read_data_arr); 
        precice_bcs[i]->update_flag = true;
    }
}

void PreciceAdapter::WriteData(const mfem::Vector T, const ConductivityModel* cond_model)
{
    for (int  i = 0; i < num_bcs; i++)
    {
        precice_bcs[i]->RetrieveWriteData(T, cond_model);
        interface->writeBlockScalarData(precice_bcs[i]->write_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->write_data_arr);
    }
}

void PreciceAdapter::SaveOldState(const Vector T)
{
    old_state_T = T;
}
void PreciceAdapter::ReloadOldState(Vector& T) const
{
    T = old_state_T;
}

PreciceAdapter::~PreciceAdapter()
{
    delete[] precice_bcs;
}