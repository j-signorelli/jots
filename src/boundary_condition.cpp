#include "boundary_condition.hpp"

using namespace std;
using namespace mfem;
using namespace precice;

string UniformConstantIsothermalBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Isothermal --- Value = " << uniform_value;
    return sstm.str();
}

string UniformConstantHeatFluxBC::GetInitString() const
{   
    stringstream sstm;
    sstm << "Heat Flux --- Value = " << uniform_value;
    return sstm.str();
}

string UniformSinusoidalIsothermalBC::GetInitString() const
{
    stringstream sstm;
    sstm << "Sinusoidal Isothermal --- T = " << amplitude << "*sin(" << ang_freq << "t + " << phase << ") + " << vert_shift;
    return sstm.str();
}

string UniformSinusoidalHeatFluxBC::GetInitString() const
{
    stringstream sstm;
    sstm << "Sinusoidal Heat Flux --- q_wall = " << amplitude << "*sin(" << ang_freq << "t + " << phase << ") + " << vert_shift;
    return sstm.str();
}

string PreciceIsothermalBC::GetInitString() const
{
    stringstream sstm;
    sstm << "preCICE Isothermal --- Mesh: " << mesh_name << " --- Default Value = " << default_value;
    return sstm.str();
}

string PreciceHeatFluxBC::GetInitString() const
{
    stringstream sstm;
    sstm << "preCICE Heat Flux --- Mesh: " << mesh_name << " --- Default Value = " << default_value;
    return sstm.str();
}

UniformConstantBC::UniformConstantBC(const int attr, const double in_value)
: BoundaryCondition(attr),
  uniform_value(in_value)
{
    // Initialize coefficient
    coeff = new ConstantCoefficient(uniform_value);
}
        
UniformSinusoidalBC::UniformSinusoidalBC(const int attr, const double in_amp, const double in_angfreq, const double in_phase, const double in_vert)
: BoundaryCondition(attr),
  amplitude(in_amp),
  ang_freq(in_angfreq),
  phase(in_phase),
  vert_shift(in_vert)
{   
    // Define the function
    function<double(const Vector&, double)> TDF = [=](const Vector&x, double t) -> double { return amplitude*sin(ang_freq*t + phase) + vert_shift;};
    
    // Initialize the coefficient
    coeff = new FunctionCoefficient(TDF);

}

void UniformSinusoidalBC::UpdateCoeff(const double time)
{ 
    coeff->SetTime(time);
}

PreciceBC::PreciceBC(const int attr, ParFiniteElementSpace& f, const string in_mesh, const double in_value, const string in_read, const string in_write) 
: BoundaryCondition(attr),
  mesh_name(in_mesh),
  default_value(in_value),
  read_data_name(in_read),
  write_data_name(in_write),
  dim(f.GetMesh()->Dimension()),
  update_flag(false),
  coords(nullptr),
  read_data_arr(nullptr),
  write_data_arr(nullptr),
  vertex_ids(nullptr),
  coeff_gf(&f)

{
    // Get and save coordinate array of all vertices, all bdr elements, and all bdr dof indices
    const FiniteElement *fe;
    ElementTransformation *transf;

    Array<double> coords_temp(0);

    // Loop through all boundary elements
    for (size_t i = 0; i < f.GetNBE(); i++)
    {
        // Skip over elements not on this BCs boundary
        if (f.GetBdrAttribute(i) != attr)
            continue;

        // Save boundary element index
        bdr_elem_indices.Append(i);
        fe = f.GetBE(i);
        transf = f.GetBdrElementTransformation(i);
        const IntegrationRule &ir = fe->GetNodes(); // Get the nodal DOFs in reference coords

        Array<int> bdr_elem_dofs;
        f.GetBdrElementDofs(i, bdr_elem_dofs);

        // Append bdr_elem_dofs to bdr_dof_indices member variable
        bdr_dof_indices.Append(bdr_elem_dofs);

        // Loop through all bdr element dofs
        // In this case: fe->GetDof == ir.GetNPoints
        for (size_t j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            // Set x,y,z of each dof into respective arrays
            Vector coord(3);
            transf->Transform(ip, coord); // Transform the integration point from reference to physical coords
            coords_temp.Append(coord[0]);
            coords_temp.Append(coord[1]);
            if (dim > 2)
            {
                coords_temp.Append(coord[2]);
            }
        }
    }

    // Get coordinates as single double* array
    coords = new double[coords_temp.Size()];
    for (size_t i = 0; i < coords_temp.Size(); i++)
        coords[i] = coords_temp[i];

    num_dofs = bdr_dof_indices.Size();

    // Create arrays for vertex ids, read data, and write data
    vertex_ids = new int[num_dofs];
    read_data_arr = new double[num_dofs];
    write_data_arr = new double[num_dofs];

    // Initialize coefficient
    coeff = new GridFunctionCoefficient(&coeff_gf);
}

void PreciceBC::UpdateCoeff(const double time)
{   
    if (update_flag) // If flagged for update by adapter, update gf linked to coefficient
    {
        // Update the GridFunction bdr values
        coeff_gf.SetSubVector(bdr_dof_indices, read_data_arr);
        update_flag = false;
    }
}

void PreciceIsothermalBC::RetrieveWriteData(const ParGridFunction &u_gf)
{
    // Need to get gradient of solution --> more involved!
    const FiniteElement *fe;
    ElementTransformation *transf;
    const ParFiniteElementSpace& fespace = *(u_gf.ParFESpace());

    // Loop through boundary elements
    int node_index = 0;
    for (size_t i = 0; i < bdr_elem_indices.Size(); i++)
    {
        fe = fespace.GetBE(bdr_elem_indices[i]);
        transf = fespace.GetBdrElementTransformation(bdr_elem_indices[i]);

        const IntegrationRule &ir = fe->GetNodes(); // Get nodes of the bdr element

        // Loop through this bdr element's DOFs
        // Again, fe->GetDof == ir.GetNPoints
        for (size_t j = 0; j < fe->GetDof(); j++)
        {
            const IntegrationPoint &ip = ir.IntPoint(j);
            transf->SetIntPoint(&ip); // Set integration pt to get appropriate Jacobian

            // Get the wall normal vector
            Vector normal(transf->Jacobian().Height());
            CalcOrtho(transf->Jacobian(), normal);

            // Get the gradient of temperature at the given nodal value for this bdr element
            Vector grad_T(normal.Size());
            u_gf.GetGradient(*transf, grad_T);

            // Get local thermal conductivity (NOTE: assumed here again of isotropic thermal conductivity)
            double k_loc = k_prop.GetLocalValue(u_gf.GetValue(*transf, ip));

            // Calculate + set value of heat
            write_data_arr[node_index] = - k_loc * (grad_T * normal) / normal.Norml2();
            node_index++;
        }
    }
}


void PreciceHeatFluxBC::RetrieveWriteData(const ParGridFunction &u_gf)
{
    // Just get solution data @ bdr
    u_gf.GetSubVector(bdr_dof_indices, write_data_arr);
}

PreciceBC::~PreciceBC()
{
    delete[] coords;
    delete[] vertex_ids;
    delete[] read_data_arr;
    delete[] write_data_arr;
}
