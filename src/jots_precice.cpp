#include "jots_precice.hpp"

using namespace std;
using namespace precice;
using namespace mfem;

const string PreciceAdapter::cowid = precice::constants::actionWriteInitialData();
const string PreciceAdapter::cowic = precice::constants::actionWriteIterationCheckpoint();
const string PreciceAdapter::coric = precice::constants::actionReadIterationCheckpoint();

PreciceAdapter::PreciceAdapter(string in_part_name, string in_config, const int r, const int s, MPI_Comm comm)
: interface(in_part_name, in_config, r, s, &comm),
  precice_bcs(nullptr),
  dim(interface.getDimensions())
{

}

void PreciceAdapter::SetPreciceBCs(BoundaryCondition** in_bcs, vector<int> precice_bc_indices)
{
    num_bcs = precice_bc_indices.size();
    precice_bcs = new PreciceBC*[num_bcs];

    // Add BCs to list
    for (size_t i = 0; i < num_bcs; i++)
    {
        precice_bcs[i] = dynamic_cast<PreciceBC*>(in_bcs[precice_bc_indices[i]]);

        // Set mesh ID
        int mesh_id = interface.getMeshID(precice_bcs[i]->mesh_name);
        precice_bcs[i]->mesh_id = mesh_id;

        // Set mesh vertices
        interface.setMeshVertices(mesh_id, precice_bcs[i]->num_dofs, precice_bcs[i]->coords, precice_bcs[i]->vertex_ids);

        // Set read/write data IDs
        precice_bcs[i]->read_data_id = interface.getDataID(precice_bcs[i]->read_data_name, mesh_id);
        precice_bcs[i]->write_data_id = interface.getDataID(precice_bcs[i]->write_data_name, mesh_id);
    }

}


void PreciceAdapter::GetReadData()
{
    for (size_t i = 0; i < num_bcs; i++)
    {
        interface.readBlockScalarData(precice_bcs[i]->read_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->read_data_arr); 
        
        // If Neumann, negate incoming values
        if (!precice_bcs[i]->IsEssential())
        {
            for (int j = 0; j < precice_bcs[i]->num_dofs; j++)
                precice_bcs[i]->read_data_arr[j] *= -1;
        }
        
        precice_bcs[i]->update_flag = true;
    }
}

void PreciceAdapter::WriteData(const ParGridFunction &u_gf)
{
    for (size_t i = 0; i < num_bcs; i++)
    {
        precice_bcs[i]->RetrieveWriteData(u_gf);
        interface.writeBlockScalarData(precice_bcs[i]->write_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->write_data_arr);
    }
}

void PreciceAdapter::SaveOldState(const Vector u)
{
    old_state_u = u;
}
void PreciceAdapter::ReloadOldState(Vector& u) const
{
    u = old_state_u;
}

PreciceAdapter::~PreciceAdapter()
{
    delete[] precice_bcs;
}