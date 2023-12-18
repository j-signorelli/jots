#include "jots_nlfs.hpp"


NonlinearJOTSDiffusionIntegrator::NonlinearJOTSDiffusionIntegrator(MaterialProperty& k_, ParFiniteElementSpace* fespace_)
: k(k_),
  fespace(*fespace_),
  u_gf(fespace_),
  u_coeff(&u_gf),
  dkdu_times_u(k.GetDCoeffRef(), u_coeff),
  diff(k.GetCoeffRef()),
  diff_d(dkdu_times_u)
{

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Update the coefficient k
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);

    // Very simply can just use DiffusionIntegrator::AssembleElementVector,
    // as the coefficient is re-evaluated every call
    diff.AssembleElementVector(el, Tr, elfun, elvect);

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    // Update both coefficients as both required to be updated
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);
    k.UpdateDCoeff(elfun, dofs);
    
    // Update u GridFunction associated with dkdu * u
    u_gf.SetSubVector(dofs, elfun);

    // Sum of usual diffusion term plus nonlinear term
    // Very simply can just use DiffusionIntegrator::AssembleElementMatrix,
    // as the coefficient is re-evaluated every call
    diff.AssembleElementMatrix(el, Tr, elmat);

    // Can use DiffusionIntegrator but with lambda = k'(u)u
    DenseMatrix elmat_d;
    diff_d.AssembleElementMatrix(el, Tr, elmat_d);

    elmat += elmat_d;


}