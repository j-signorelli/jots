#include "jots_precice.hpp"

using namespace std;
using namespace precice;
using namespace mfem;

const string JOTSSolverInterface::cowid = precice::constants::actionWriteInitialData();
const string JOTSSolverInterface::cowic = precice::constants::actionWriteIterationCheckpoint();
const string JOTSSolverInterface::coric = precice::constants::actionReadIterationCheckpoint();

JOTSSolverInterface::JOTSSolverInterface(string in_part_name, string in_config, const int r, const int s, MPI_Comm comm)
: SolverInterface(in_part_name, in_config, r, s, &comm),
  precice_bcs(nullptr)
{

}

void JOTSSolverInterface::SetPreciceBCs(BoundaryCondition** in_bcs, vector<int> precice_bc_indices)
{
    num_bcs = precice_bc_indices.size();
    precice_bcs = new PreciceBC*[num_bcs];

    // Add BCs to list
    for (size_t i = 0; i < num_bcs; i++)
    {
        precice_bcs[i] = dynamic_cast<PreciceBC*>(in_bcs[precice_bc_indices[i]]);

        // Set mesh ID
        int mesh_id = getMeshID(precice_bcs[i]->mesh_name);
        precice_bcs[i]->mesh_id = mesh_id;

        // Set mesh vertices
        setMeshVertices(mesh_id, precice_bcs[i]->num_dofs, precice_bcs[i]->coords, precice_bcs[i]->vertex_ids);

        // Set read/write data IDs
        precice_bcs[i]->read_data_id = getDataID(precice_bcs[i]->read_data_name, mesh_id);
        precice_bcs[i]->write_data_id = getDataID(precice_bcs[i]->write_data_name, mesh_id);
    }

}


void JOTSSolverInterface::GetReadData()
{
    for (size_t i = 0; i < num_bcs; i++)
    {
        readBlockScalarData(precice_bcs[i]->read_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->read_data_arr); 
        
        // If Neumann, negate incoming values
        if (!precice_bcs[i]->IsEssential())
        {
            for (int j = 0; j < precice_bcs[i]->num_dofs; j++)
                precice_bcs[i]->read_data_arr[j] *= -1;
        }
        
        precice_bcs[i]->update_flag = true;
    }
}

void JOTSSolverInterface::WriteData(const ParGridFunction &u_gf)
{
    for (size_t i = 0; i < num_bcs; i++)
    {
        precice_bcs[i]->RetrieveWriteData(u_gf);
        writeBlockScalarData(precice_bcs[i]->write_data_id, precice_bcs[i]->num_dofs, precice_bcs[i]->vertex_ids, precice_bcs[i]->write_data_arr);
    }
}

void JOTSSolverInterface::SaveOldState(const Vector u)
{
    old_state_u = u;
}
void JOTSSolverInterface::ReloadOldState(Vector& u) const
{
    u = old_state_u;
}

JOTSSolverInterface::~JOTSSolverInterface()
{
    delete[] precice_bcs;
}