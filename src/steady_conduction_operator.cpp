#include "steady_conduction_operator.hpp"


using namespace mfem;

NonlinearJOTSDiffusionIntegrator(const MaterialProperty& k_, ParFiniteElementSpace* fespace_)
: k(k_),
  fespace(fespace_),
  f(fespace),
  df(fespace),
  f_coeff(*f),
  df_coeff(*df),
  integ_vec(f_coeff),
  integ_grad(df_coeff)
{

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Loop through DOFs, update elvec w/ MaterialProperty::GetLocalValue results
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    // Write code to get vector f(elfun)
    f.SetSubVector(dofs, __INSERT_f(elfun)__HERE);
    integ_vec.AssembleRHSElementVect(el, Tr, elvect);
}

void NonlinearJOTSDiffusionIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    //**Need to understand what exactly this is.....
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    // Write code to get vector df(elfun)
    df.SetSubVector(dofs, __INSERT_df(elfun)__HERE);
    integ_grad.AssembleElementMatrix(el, Tr, elmat);
}

SteadyConductionOperator::SteadyConductionOperator()
: 
{

}