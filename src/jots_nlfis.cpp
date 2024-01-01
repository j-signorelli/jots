#include "jots_nlfis.hpp"


JOTSNonlinearDiffusionIntegrator::JOTSNonlinearDiffusionIntegrator(ParFiniteElementSpace* fespace_, MaterialProperty& k_)
: fespace(*fespace_),
  lambda(lambda_),
  dlambdadu(dlambdadu_)
  u_gf(fespace_),
  grad_u_coeff(&u_gf),
  dlambdadu_times_grad_u(dlambdadu_, grad_u_coeff),
  term1(dlambdadu_times_grad_u),
  term2(lambda)
{

}

void JOTSNonlinearDiffusionIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Very simply can just use DiffusionIntegrator::AssembleElementVector,
    // as the coefficient is re-evaluated every call
    term2.AssembleElementVector(el, Tr, elfun, elvect);

}

void JOTSNonlinearDiffusionIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    // Update u GridFunction associated with dkdu * grad u
    fespace.GetElementDofs(Tr.ElementNo, dofs);
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


NonlinearJOTSMassIntegrator::NonlinearJOTSMassIntegrator(ParFiniteElementSpace* fespace_, MaterialProperty& rho_, MaterialProperty& C_)
: fespace(*fespace_),
  rho(rho_),
  C(C_),
  u_gf(fespace_),
  u_coeff(&u_gf),
  drhodu_C(rho.GetDCoeffRef(), C.GetCoeffRef()),
  rho_dCdu(rho.GetCoeffRef(), C.GetDCoeffRef()),
  mat_prop_coeff(drhodu_C, rho_dCdu),
  mat_prop_coeff_times_u(mat_prop_coeff, u_coeff),
  term1(mat_prop_coeff_times_u),
  term2(rho_C)
{

}

NonlinearJOTSMassIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    // Update both coefficients
    rho.UpdateCoeff(elfun, dofs);
    C.UpdateCoeff(elfun, dofs);

    term2.AssembleElementVector(el, Tr, elfun, elvect);
}

NonlinearJOTSMassIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    // Update both coefficients + their derivatives
    rho.UpdateCoeff(elfun, dofs);
    C.UpdateCoeff(elfun, dofs);
    rho.UpdateDCoeff(elfun, dofs);
    C.UpdateDCoeff(elfun, dofs);


    // Update u GridFunction associated with mat_prop_coeff* u
    u_gf.SetSubVector(dofs, elfun);


    term2.AssembleElementMatrix(el, Tr, elmat);
    DenseMatrix elmat_d;
    term1.AssembleElementMatrix(el, Tr, elmat_d);
    elmat += elmat_d;
}