#include "boundary_condition.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

string UniformIsothermalBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Isothermal --- Value: " << uniform_value;
    return sstm.str();
}

string UniformHeatFluxBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Heat Flux --- Value: " << uniform_value;
    return sstm.str();
}


preCICEBC::preCICEBC(int attr, BOUNDARY_CONDITION in_type, SolverInterface* in, ParGridFunction* in_T_gf, bool is_restart, string mesh_name, double initial_value); 
: BoundaryCondition(attr, in_type), 
  fespace(*in_T_gf->FESpace()),
  interface(in),
  T_gf(in_T_gf),
  restart(is_restart),
  coeff_gf(f),
  boundary_dofs(0),
  coeff_values(0)
{
    // Method from GridFunction::AccumulateAndCountBdrValues used

    const FiniteElement *fe;
    ElementTransformation *transf;

    int dim = fespace->GetMesh()->Dimension();
    int num_dofs;

    Array<double> coords_temp(0);

    // Loop through all boundary elements
    for (int i = 0; i < fespace->GetNBE(); i++)
    {   
        // Skip over elements not on this BCs boundary
        if (fespace->GetBdrAttribute(i) != 1)
            continue;

        fe = fespace->GetBE(i); // Get bdr element
        num_dofs = fe->GetDof();

        transf = fespace->GetBdrElementTransformation(i);
        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the element

        Array<int> bdr_elem_dofs;
        fespace->GetBdrElementDofs(i, bdr_elem_dofs); // Get bdr element dofs

        // Append bdr_elem_dofs to boundary_dofs member variable
        boundary_dofs.Append(bdr_elem_dofs);

        // Loop through all bdr element dofs
        for (int j = 0; j < num_dofs; j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip);
            
            // Set x,y,z of each dof into respective arrays
            Vector coord(3);
            transf->Transform(ip, coord);
            coords_temp.Append(coord[0]);
            coords_temp.Append(coord[1]);
            if (dim > 2)
            {
                coords_temp.Append(coord[2]);
            }

        }
    }

    /* TODO: Remove if unnecessary later. (If ghost nodes values are correct, then below unnecessary)
    // Now get all bdr tdofs and create new array/vector of only true dofs from above
    Array<int> bdr_tdofs;
    fespace->GetBoundaryTrueDofs(bdr_tdofs, 1);
    bool tdof = false;
    for (int i = 0; i < boundary_dofs.Size(); i++)
    {
        for (int j = 0; j < bdr_tdofs.Size(); j++)
        {
            if (boundary_dofs[i] == bdr_tdofs[j])
            {
                tdof = true;
                j = bdr_tdofs.Size();
            }
        }
        if (!tdof)
            cout << "NON-TRUE DOF: " << boundary_dofs[i] << endl;
    }
    */

    // Get coordinates as single double* array
    coords = new double[coords_mfem.Size()];
    for (int i = 0; i < coords_fem.Size(); i++)
        coords[i] = coords_mfem[i];

    // Get mesh ID
    meshID = interface->getMeshID(mesh_name);

    // Set mesh vertices
    interface->SetMeshVertices(meshID, boundary_dofs.Size(), coords, vertexIDs);

    // Get read + write data IDs
    readDataID = interface->getDataID(GetReadDataName(), meshID);
    writeDataID = interface->getDataID(GetWriteDataName(), meshID);

    // Set coeff_values to initialization value for now
    // If data sent from other participant, this will be updated in first call to UpdateCoeff
    coeff_values.SetSize(boundary_dofs.Size());
    coeff_values = initial_value;
}

double* preCICEBC::GetTemperatures(ParGridFunction* in_T_gf, Array<int> bdr_dofs)
{

}

double* preCICEBC::GetWallHeatFlux(ParGridFunction* in_T_gf, Array<int> bdr_dofs)
{

}

//void preCICEBC::InitCoefficient()
//{
    // This is where restart/nonrestart becomes important.
    // If restart, use current GF values
    // If nonrestart, need to apply BC or something first
    // Depends on coefficient

//}
