#include "jots_iterator.hpp"

using namespace mfem;

JOTSIterator::JOTSIterator(ParFiniteElementSpace& f_, const Config &in_config, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, int bc_count)
: fespace(f_),
  b(&f_),
  b_vec(f_.GetTrueVSize()),
  neumann_coeff(fespace.GetVDim()),
  lin_solver(nullptr)
{
    int dim = fespace.GetVDim();

    // Set linear solver + preconditioner
    lin_solver = Factory::GetSolver(in_config.GetSolverLabel(), fespace.GetComm());
    lin_solver->iterative_mode = false;
    lin_prec.SetType(Factory::GetPrec(in_config.GetPrecLabel()));
    lin_solver->SetPreconditioner(lin_prec);
    lin_solver->SetAbsTol(in_config.GetAbsTol());
    lin_solver->SetRelTol(in_config.GetRelTol());
    lin_solver->SetMaxIter(in_config.GetMaxIter());
    lin_solver->SetPrintLevel(Factory::CreatePrintLevel(in_config.GetLinSolPrintLevel()));


    // Get the essential DOFs component-wise
    // Set Neumann coefficient component-wise
    for (int comp = 0; comp < dim; comp++)
    {
        PWCoefficient* comp_neumann_coeff = new PWCoefficient();

        // Set the list of Dirichlet (essential) DOFs (for this component)
        Array<int> dbc_bdr_comp(bc_count);
        dbc_bdr_comp = 0; // Start w/ all attributes set to non-essential = 0

        for (int i = 0; i < bc_count; i++)
        {
            if (in_bcs[i*dim+comp]->IsEssential())
                dbc_bdr_comp[i] = 1;
            else
            {
                comp_neumann_coeff->UpdateCoefficient(in_bcs[i*dim+comp]->GetBdrAttr(), in_bcs[i*dim+comp]->GetCoeffRef());
            }
        }
        // Get the essential true dofs, given the Dirichlet boundaries FOR THIS COMPONENT
        Array<int> ess_tdof_list_comp;
        fespace.GetEssentialTrueDofs(dbc_bdr_comp, ess_tdof_list_comp, fespace.GetVDim() == 1 ? -1 : comp);
        // Add to full list
        ess_tdof_list.Append(ess_tdof_list_comp);

        // Add this component Neumann coefficient
        neumann_coeff.Set(comp, comp_neumann_coeff, true); // VectorArrayCoefficient will take ownership of it + delete it

    }

    // Add Neumann term
    b.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(neumann_coeff));
    
    // Initialize Neumann LF and vector
    UpdateNeumann();
}

void JOTSIterator::UpdateNeumann()
{
    b.Assemble();
    b.ParallelAssemble(b_vec);
}