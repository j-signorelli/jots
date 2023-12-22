#include "jots_nlfs.hpp"


NonlinearJOTSDiffusionIntegrator::NonlinearJOTSDiffusionIntegrator(MaterialProperty& k_, ParFiniteElementSpace* fespace_)
: k(k_),
  fespace(*fespace_),
  u_gf(fespace_),
  grad_u_coeff(&u_gf),
  dkdu_times_grad_u(k.GetDCoeffRef(), grad_u_coeff),
  term1(dkdu_times_grad_u),
  term2(k.GetCoeffRef())
{

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Update the coefficient k
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);

    // Very simply can just use DiffusionIntegrator::AssembleElementVector,
    // as the coefficient is re-evaluated every call
    term2.AssembleElementVector(el, Tr, elfun, elvect);

}

void NonlinearJOTSDiffusionIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    // Update both coefficients as both required to be updated
    fespace.GetElementDofs(Tr.ElementNo, dofs);
    k.UpdateCoeff(elfun, dofs);
    k.UpdateDCoeff(elfun, dofs);
    
    // Update u GridFunction associated with dkdu * grad u
    u_gf.SetSubVector(dofs, elfun);

    // Sum of usual diffusion term plus nonlinear term
    // Very simply can just use DiffusionIntegrator::AssembleElementMatrix,
    // as the coefficient is re-evaluated every call
    term2.AssembleElementMatrix(el, Tr, elmat);

    DenseMatrix elmat_d;
    term1.AssembleElementMatrix(el, Tr, elmat_d);

    elmat_d.Neg();
    elmat += elmat_d;


}

/*
NonlinearJOTSMassIntegrator::NonlinearJOTSMassIntegrator(MaterialProperty& rho_, MaterialProperty& C_, ParFiniteElementSpace* fespace_)
: rho(rho_),
  C(C_),
  fespace(*fespace_),
  u_gf(fespace_),
  u_coeff(&u_gf),
  rho_C(rho.GetLocalValue(0), C.GetCoeffRef()),
  rho_dCdu(rho.GetLocalValue(0), C.GetDCoeffRef()),
  rho_dCdu_times_u(rho_dCdu, u_coeff),
  mass(rho_C),
  mass_d(rho_dCdu_times_u)
{

}
*/