#include "boundary_condition.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

string UniformConstantIsothermalBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Isothermal --- Value: " << uniform_value;
    return sstm.str();
}

string UniformConstantHeatFluxBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Heat Flux --- Value: " << uniform_value;
    return sstm.str();
}

string UniformSinusoidalIsothermalBC::GetInitString() const
{
    stringstream sstm;
    sstm << "Sinusoidal Isothermal --- T = " << amplitude << "*sin(" << ang_freq << "t + " << phase << ") + " << vert_shift;
    return sstm.str();
}

string PreciceIsothermalBC::GetInitString() const
{
    stringstream sstm;
    sstm << "preCICE Isothermal --- Mesh: " << mesh_name << " --- Default Value: " << default_value;
    return sstm.str();
}

string PreciceHeatFluxBC::GetInitString() const
{
    stringstream sstm;
    sstm << "preCICE Heat Flux --- Mesh: " << mesh_name << " --- Default Value: " << default_value;
    return sstm.str();
}
void UniformSinusoidalIsothermalBC::InitCoefficient()
{   
    function<double(const Vector&, double)> TDF = [=](const Vector&x, double t) -> double { return amplitude*sin(ang_freq*t + phase) + vert_shift;};
    coeff = new FunctionCoefficient(TDF);//
}
void UniformSinusoidalIsothermalBC::UpdateCoeff()
{ 
    coeff->SetTime(time_ref);
}

PreciceBC::PreciceBC(const int attr, const BOUNDARY_CONDITION in_type, ParFiniteElementSpace& f, const string in_mesh, const bool is_restart, const double in_value, const string in_read, const string in_write) 
: BoundaryCondition(attr, in_type),
  fespace(f),
  mesh_name(in_mesh),
  restart(is_restart),
  default_value(in_value),
  read_data_name(in_read),
  write_data_name(in_write),
  dim(f.GetMesh()->Dimension()),
  coords(nullptr),
  vertex_ids(nullptr),
  read_data_arr(nullptr),
  write_data_arr(nullptr),
  update_flag(false),
  bdr_elem_indices(0),
  bdr_dof_indices(0)
  //coeff_values(0),

{
    // Method from GridFunction::AccumulateAndCountBdrValues used
    // Get and save coordinate array of all vertices, all bdr elements, and all bdr dof indices

    //Get all bdr tdofs and create new array/vector of only true dofs from above
    //Array<int> bdr_tdofs;
    //fespace->GetBoundaryTrueDofs(bdr_tdofs, 1);

    const FiniteElement *fe;
    ElementTransformation *transf;

    Array<double> coords_temp(0);

    // Loop through all boundary elements
    for (int i = 0; i < fespace.GetNBE(); i++)
    {   
        // Skip over elements not on this BCs boundary
        if (fespace.GetBdrAttribute(i) != attr)
            continue;
        
        // Save boundary element index
        bdr_elem_indices.Append(i);

        fe = fespace.GetBE(i); // Get bdr element used at this boundary face i
        transf = fespace.GetBdrElementTransformation(i); // Get bdr elem transformation - the actual definition of this particular BE

        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the bdr element

        Array<int> bdr_elem_dofs;
        fespace.GetBdrElementDofs(i, bdr_elem_dofs); // Get bdr element dofs

        // Append bdr_elem_dofs to boundary_dofs member variable
        bdr_dof_indices.Append(bdr_elem_dofs);

        // Loop through all bdr element dofs
        for (int j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip); // TODO: is this line unnecessary? See docs

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
    /*
    //Now get all bdr tdofs and create new array/vector of only true dofs from above
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
        //if (!tdof)
        //    cout << "NON-TRUE DOF: " << boundary_dofs[i] << endl;
    }
    */
   // ^ TODO After functioning in serial

    // Get coordinates as single double* array
    coords = new double[coords_temp.Size()];
    for (int i = 0; i < coords_temp.Size(); i++)
        coords[i] = coords_temp[i];

    num_dofs = bdr_dof_indices.Size();

    // Create arrays for everything else needed
    vertex_ids = new int[num_dofs];
    read_data_arr = new double[num_dofs];
    write_data_arr = new double[num_dofs];
    

    // Set coeff_values to initialization value for now
    // If data sent from other participant, this will be updated in first call to UpdateCoeff
    coeff_dof_values.SetSize(num_dofs);
    coeff_dof_values = default_value;

    // Create grid function for coefficient
    coeff_gf = new ParGridFunction(&fespace);

    // Create auxiliary GF on fespace
    temp_gf = new ParGridFunction(&fespace);
}

void PreciceBC::InitCoefficient()
{   
    coeff = new GridFunctionCoefficient(coeff_gf);
    // Do not need to reset coeff_gf to coeff, just update it!
    // (coeff takes ptr to coeff_gf) 
}

void PreciceBC::UpdateCoeff()
{   

    if (update_flag) // If flagged for update by adapter, update gf linked to coefficient
    {
        // Update the GridFunction bdr values
        coeff_gf->SetSubVector(bdr_dof_indices, read_data_arr);
        update_flag = false;
    }
}

// TODO: Can I somehow not send cond_model pointlessly? Maybe check if essential or something??
void PreciceIsothermalBC::RetrieveInitialWriteData(const mfem::Vector T, const ConductivityModel* cond_model)
{
    // For isothermal wall, sending heat flux
    //      If not restart, must project coeff onto initialization first, then get heat flux, then send
    //      If restart, do nothing
    temp_gf->SetFromTrueDofs(T);
    if (!restart) // TODO: I don't think this is necessary. See what SU2 does. it may not enforce BCs for initialization.
    {
        ConstantCoefficient temp_coeff(default_value);
        temp_gf->ProjectCoefficient(temp_coeff, bdr_dof_indices);// TODO: can I just setsubvector for H1??
    }
    GetBdrWallHeatFlux(temp_gf, cond_model, bdr_elem_indices, write_data_arr);
}

void PreciceIsothermalBC::RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model)
{
    temp_gf->SetFromTrueDofs(T);
    GetBdrWallHeatFlux(temp_gf, cond_model, bdr_elem_indices, write_data_arr);
}

void PreciceHeatFluxBC::RetrieveInitialWriteData(const mfem::Vector T, const ConductivityModel* cond_model)
{
    // For heatflux wall, sending temperature
    //      If restart, get temperature from state and send
    //      If not restart, not at all any difference.
    RetrieveWriteData(T, cond_model);
}

void PreciceHeatFluxBC::RetrieveWriteData(const mfem::Vector T, const ConductivityModel* cond_model)
{
    temp_gf->SetFromTrueDofs(T);
    GetBdrTemperatures(temp_gf, bdr_elem_indices, write_data_arr);
}

void PreciceBC::GetBdrTemperatures(const ParGridFunction* T_gf, const Array<int> in_bdr_elem_indices, double* nodal_temperatures)
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
            // TODO: Can I just do GetSubVector on T for H1?
        }
    }
}

void PreciceBC::GetBdrWallHeatFlux(const mfem::ParGridFunction* T_gf, const ConductivityModel* in_cond, const mfem::Array<int> in_bdr_elem_indices, double* nodal_wall_heatfluxes)
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
            double k = in_cond->GetLocalConductivity(T_gf->GetValue(*transf, ip));

            // Calculate + set value of heat
            nodal_wall_heatfluxes[nodal_index] = k * (grad_T * normal) / normal.Norml2();       
            nodal_index++;
        }
    }
}

PreciceBC::~PreciceBC()
{
    delete[] coords;
    delete[] vertex_ids;
    delete[] read_data_arr;
    delete[] write_data_arr;
    delete coeff_gf;
    delete temp_gf;
}
