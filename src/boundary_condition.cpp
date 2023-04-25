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


preCICEBC::preCICEBC(int attr, BOUNDARY_CONDITION in_type, SolverInterface* in, const ParGridFunction* in_T_gf, ConductivityModel* in_cond, bool is_restart, string mesh_name, double initial_value, string read_data_name, string write_data_name) 
: BoundaryCondition(attr, in_type),
  interface(in),
  fespace(in_T_gf->ParFESpace()),
  T_gf(in_T_gf),
  cond_model(in_cond),
  readDataName(read_data_name),
  writeDataName(write_data_name),
  restart(is_restart),
  bdr_elem_indices(0),
  bdr_dof_indices(0),
  coeff_gf(in_T_gf->ParFESpace()),
  coeff_values(0)
{
    // Method from GridFunction::AccumulateAndCountBdrValues used

    const FiniteElement *fe;
    ElementTransformation *transf;

    int dim = fespace->GetMesh()->Dimension();

    Array<double> coords_temp(0);

    // Loop through all boundary elements
    for (int i = 0; i < fespace->GetNBE(); i++)
    {   
        // Skip over elements not on this BCs boundary
        if (fespace->GetBdrAttribute(i) != attr)
            continue;
        
        // Save boundary element index
        bdr_elem_indices.Append(i);

        fe = fespace->GetBE(i); // Get bdr element used at this boundary face i
        transf = fespace->GetBdrElementTransformation(i); // Get bdr elem transformation - the actual definition of this particular BE

        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the bdr element

        Array<int> bdr_elem_dofs;
        fespace->GetBdrElementDofs(i, bdr_elem_dofs); // Get bdr element dofs

        // Append bdr_elem_dofs to boundary_dofs member variable
        bdr_dof_indices.Append(bdr_elem_dofs);

        // Loop through all bdr element dofs
        for (int j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip); // TODO: Why is this necessary for Transform?? Should see if relevant for me - See Coefficient::Eval documentation
            
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
    coords = new double[coords_temp.Size()];
    for (int i = 0; i < coords_temp.Size(); i++)
        coords[i] = coords_temp[i];

    // Get mesh ID
    meshID = interface->getMeshID(mesh_name);

    // Set mesh vertices
    interface->setMeshVertices(meshID, bdr_dof_indices.Size(), coords, vertexIDs);

    // Get read + write data IDs
    readDataID = interface->getDataID(readDataName, meshID);
    writeDataID = interface->getDataID(writeDataName, meshID);

    // Create arrays for read/write of preCICE data
    readDataArr = new double[bdr_dof_indices.Size()];
    writeDataArr = new double[bdr_dof_indices.Size()];
    // Set coeff_values to initialization value for now
    // If data sent from other participant, this will be updated in first call to UpdateCoeff
    coeff_values.SetSize(bdr_dof_indices.Size());
    coeff_values = initial_value;

    /*
    // DEBUG CHECK:
    cout << "Printing nodal values as debug check: " << endl;
    preCICEBC::SetBdrTemperatures(T_gf, bdr_elem_indices, readDataArr);
    for (int i = 0; i < bdr_dof_indices.Size(); i++)
        cout << "BDR DOF " << bdr_dof_indices[i] << ": " << readDataArr[i] << " K\n";
    */
}

void preCICEBC::GetBdrTemperatures(const ParGridFunction* T_gf, const Array<int> in_bdr_elem_indices, double* nodal_temperatures)
{   
    ParFiniteElementSpace* fespace = T_gf->ParFESpace();

    const FiniteElement *fe;
    ElementTransformation *transf;

    int nodal_index = 0;
    for (int i = 0; i < in_bdr_elem_indices.Size(); i++)
    {
        fe = fespace->GetBE(i); // Get FiniteElement
        transf = fespace->GetBdrElementTransformation(i); //Get this associated ElementTransformation from the mesh object
        // ^Above is same as just calling GetBdrElementTransformation from Mesh

        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the bdr element

        // Loop through this bdr element's DOFs
        for (int j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip); // TODO: See above- is this line unnecessary?

            // Set the value
            nodal_temperatures[nodal_index] = T_gf->GetValue(*transf, ip);       
            nodal_index++;
        }
    }
}

void preCICEBC::GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* cond_model, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes)
{

    ParFiniteElementSpace* fespace = T_gf->ParFESpace();

    const FiniteElement *fe;
    ElementTransformation *transf;

    int nodal_index = 0;
    for (int i = 0; i < in_bdr_elem_indices.Size(); i++)
    {
        fe = fespace->GetBE(i); // Get FiniteElement
        transf = fespace->GetBdrElementTransformation(i); //Get this associated ElementTransformation from the mesh object
        // ^Above is same as just calling GetBdrElementTransformation from Mesh

        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the bdr element

        // Loop through this bdr element's DOFs
        for (int j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip); // Set integration pt to get appropriate Jacobian

            // Get the wall normal vector
            Vector normal(transf->Jacobian().Height());
            CalcOrtho(transf->Jacobian(), normal);

            // Get the gradient of temperature
            Vector grad_T(normal.Size());
            T_gf->GetGradient(*transf, grad_T);

            // Get local thermal conductivity (NOTE: assumed here again of isotropic thermal conductivity)
            double k = cond_model->GetLocalConductivity(T_gf->GetValue(*transf, ip));

            // Calculate + set value of heat
            nodal_wall_heatfluxes[nodal_index] = grad_T * nor / nor.Norml2();       
            nodal_index++;
        }
    }
}

void preCICEBC::InitCoefficient()
{

}

void preCICEBC::UpdateCoeff()
{

}

 preCICEBC::~preCICEBC()
{
    delete[] coords;
    delete[] vertexIDs;
    delete[] readDataArr;
    delete[] writeDataArr;
}

//void preCICEBC::InitCoefficient()
//{
    // This is where restart/nonrestart becomes important.
    // If restart, use current GF values
    // If nonrestart, need to apply BC or something first
    // Depends on coefficient, but can probably setup w virtual fxns

//}
