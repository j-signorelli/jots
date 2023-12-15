#include "jots_iterator.hpp"

using namespace mfem;

JOTSIterator::JOTSIterator(ParFiniteElementSpace& f_, const BoundaryCondition* const* in_bcs, Array<int>* all_bdr_attr_markers, int bc_count)
: fespace(f_),
  b(&f_),
  b_vec_full(f_.GetTrueVSize())
{

   // Set the list of Dirichlet (essential) DOFs
   Array<int> dbc_bdr(bc_count);
   dbc_bdr = 0; // Start w/ all attributes set to non-essential = 0

   for (int i = 0; i < bc_count; i++)
   {
      // Add Neumann BCs to linear form b
      if (!in_bcs[i]->IsEssential())
      {
         b.AddBoundaryIntegrator(new BoundaryLFIntegrator(in_bcs[i]->GetCoeffRef()), all_bdr_attr_markers[i]);// Add boundary integrator to boundary
      }
      else
      {
         // Update DBC list if essential
         dbc_bdr[i] = 1;
      }
   }

   // Get the essential true dofs, given the Dirichlet boundaries. Store in ess_tdof_list
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // Initialize Neumann linear form
   UpdateNeumann();
}

void JOTSIterator::UpdateNeumann()
{
   b.Assemble();
   b.ParallelAssemble(b_vec_full);
}