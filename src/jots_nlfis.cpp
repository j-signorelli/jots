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

JOTSNonlinearNeumannIntegrator::JOTSNonlinearNeumannIntegrator(Coefficient& lambda_, Coefficient& dlambdadu_)
: lambda(lambda_),
  dlambdadu(dlambdadu_),
  vec_integ(lambda),
  grad_integ(dlambdadu)
{

}

void JOTSNonlinearNeumannIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    vec_integ.AssembleRHSElementVect(el, Tr, elvect);
}

void JOTSNonlinearNeumannIntegrator::AssembleElementGrad(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, DenseMatrix &elmat)
{
    grad_integ.AssembleElementMatrix(el, Tr, elmat);
}