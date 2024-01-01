# Unsteady Heat Transfer

## Approach #1

### Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. We assume isotropic density $\rho$, thermal conductivity $k$, and specific heat $C$ but allow for them to vary as a function of temperature: $\rho=\rho(u)$, $k=k(u)$, and $C=C(u)$
2. No heat generation

This yields:

$$ \rho C\dfrac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = 0$$

We seek a weak formulation of this allowing for nonhomogeneous Neumann and Dirichlet boundary conditions.

Given that one applies the appropriate "lifting" of the solution, we can assume homogeneous Dirichlet BCs in the derivation such that:

$$u=0 \text{ on } \partial \Omega_D$$
$$\hat{n} \cdot (k \nabla u)=g \text{ on } \partial \Omega_N \text{ given } g : \partial \Omega_N \times \mathbb{R}\rightarrow \mathbb{R}$$

where $\Omega$ is our domain and $\partial \Omega$ is the boundary of the domain. Note that $g$ takes in both a location on the boundary and is also allowed to take in a time value that is real.

So: we seek a solution that at a given point in time $u \in H^1_{\partial \Omega_D}(\Omega)$ where $H^1_{\partial \Omega_D}(\Omega) = \{u \in H^1(\Omega) : u=0 \text{ on } \partial \Omega_D \}$. Using conventional notation, the solution space is $\mathcal{V}=H^1_{\partial \Omega_D}(\Omega)$. The test space is chosen to be the same as the solution space.

Multiply the heat equation by test function $v \in H^1_{\partial \Omega_D}(\Omega)$ and integrate:

$$\int_\Omega\rho C\dfrac{\partial u}{\partial t}vd\vec{x} - \int_\Omega \nabla \cdot (k \nabla u)vd\vec{x} = 0$$

We can then apply the following Green formula to the second term:

$$\int_\Omega\nabla \cdot (k \nabla u)vd\vec{x} = -\int_\Omega (\nabla v) \cdot (k\nabla u) d\vec{x} + \int_{\partial \Omega} \hat{n} \cdot (k \nabla u)v d\vec{x}$$

Note from earlier that $ g = \hat{n} \cdot (k \nabla u)$, a specified heat flux on Neumann boundaries. We allow for this heat flux to vary in time as well as in space: $g=g(x,y,z,t)$

Now we may apply a finite element approximation to this weak formulation. This is done assuming $H^1$-conforming nodal elements. An approximate solution $u_h\in\mathcal{V}_h \sub \mathcal{V}$ is sought, where $\mathcal{V}_h$ is a finite-dimensional approximation space, which is a subset of the infinitely-dimensional solution space. Given a basis $\{\phi_j\}$ for the approximation space ($\text{span}\{\phi_j\}=\mathcal{V}_h$) and assuming that the approximate test function $v_h \in \mathcal{V}_h$, we write:
 
$$u(\vec{x},t) \approx u_h(\vec{x},t)= \sum^{N_{\mathcal{V}_h}}_j u_j(t)\phi_j(\vec{x})$$
$$v(\vec{x},t) \approx v_h(\vec{x},t)=\sum^{N_{\mathcal{V}_h}}_j v_j(t)\phi_j(\vec{x})$$


 where the dual basis $u_j$ is chosen to represent $N_{\mathcal{V}_h}$ nodal degrees of freedom in which at a given point in time: $\{u_j\}\in	\mathbb{R}^{N_{\mathcal{V}_h}}$. Going forward for simplification, all summations are presumed to be over ${N_{\mathcal{V}_h}}$ degrees of freedom: $\sum^{N_{\mathcal{V}_h}}_j=\sum_j$.

 Plugging the above in:


$$\int_\Omega \rho(u_h)C(u_h)\left(\sum_j \dfrac{d u_j(t)}{d t} \phi_j\right) \left(\sum_i v_i(t)\phi_i\right) d\vec{x} + \int_\Omega \left( k(u_h)\sum_j u_j(t)\nabla \phi_j\right)\cdot \left(\sum_i v_i(t)\nabla \phi_i\right) d\vec{x} = \int_{\partial \Omega} g\left(\sum_i v_i(t) \phi_i \right)d\vec{x}$$


Take out the summations and any non-integrated terms from the integrals:

$$\sum_j \sum_i v_i \left( \int_\Omega \rho(u_h)C(u_h)\phi_i \phi_j d\vec{x}\right)\dfrac{d u_j}{d t} + \sum_i \sum_j v_i\left( \int_\Omega k(u_h)(\nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x}\right)u_j = \sum_i v_i\left( \int_{\partial \Omega} g\phi_i d\vec{x}\right)$$

Now we may define:

$$\mathbf{M} = \left[M_{ij}\right]=\left[\int_\Omega \rho C\phi_i \phi_j d\vec{x}\right] = \text{Mass Matrix}$$

$$ \mathbf{K}= \left[K_{ij}\right] = \left[\int_\Omega k(\nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} \right] = \text{Stiffness Matrix}$$

$$\vec{N}=\left[\int_{\partial \Omega}g(\vec{x},t)\phi_i d\vec{x}\right]=\text{Neumann Term}$$

Thus the equation becomes:

$$\vec{v}^T\mathbf{M}\dfrac{d \vec{u}}{d t} + \vec{v}^T \mathbf{K}\vec{u} = \vec{v}^T \vec{N}$$

$$\vec{v}^T\left(\mathbf{M}\dfrac{d \vec{u}}{d t} + \mathbf{K}\vec{u} - \vec{N}\right)=0$$

For this inner product to be zero $\forall\vec{v}$, the vector in parentheses must be identically zero. So the equation of interest is:
$$
\mathbf{M}\dfrac{d \vec{u}}{d t} + \mathbf{K}\vec{u} - \vec{N}=\vec{0}$$

or

$$M_{ij}\dfrac{du_j}{dt} = -K_{ij}u_j + N_i$$

Note that $\vec{N}=\vec{N}(t)$. For linear cases (where all material properties are uniform) the above form with matrices holds. However, generally $\mathbf{M}=\mathbf{M}(\vec{u})$ and $\mathbf{K}=\mathbf{K}(\vec{u})$. Thus, to account for these cases, the equation must be rewritten with nonlinear operators as

$$M\left(\vec{u},\dfrac{d\vec{u}}{dt}\right) =-K(\vec{u}) + N(t)$$

where $M\left(\vec{u},\dfrac{d\vec{u}}{dt}\right)=\mathbf{M}(\vec{u})\dfrac{d\vec{u}}{dt}$ and $K(\vec{u})=\mathbf{K}(\vec{u})\vec{u}$.

### Time-Integration

JOTS makes usage of MFEM's layer-of-abstraction `TimeDependentOperator`.

#### Explicit

For explicit time integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_n$$

where $\left.\dfrac{d \vec{u}}{dt}\right|_n$ is the one that solves 

$$M\left(\vec{u}_n,\left.\dfrac{d \vec{u}}{dt}\right|_n\right) =-K(\vec{u}_n) + N(t_n)$$

Writing out the nonlinear mass operator, the following holds:

$$\left.\dfrac{d \vec{u}}{dt}\right|_n =\mathbf{M}_n^{-1}\left[-K(\vec{u}_n) + N(t_n)\right]$$

where $\mathbf{M}_n$ is the previously defined mass matrix evaluated using the solution vector at the previous timestep $t_n$.

#### Implicit

For implicit time-integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_{n+1}$$

where

$$M\left(\vec{u}_{n+1},\left.\dfrac{d \vec{u}}{dt}\right|_{n+1}\right) =-K(\vec{u}_{n+1}) + N(t_{n+1})$$

Plugging in the above equation for $\vec{u}_{n+1}$ and $t_{n+1}=t_n + \Delta t$, and writing $\vec{k}_{n+1}=\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}$

$$M(\vec{u}_n + \Delta t\vec{k}_{n+1}, \vec{k}_{n+1})=-K(\vec{u}_n + \Delta t\vec{k}_{n+1}) + N(t_n + \Delta t)$$

Rearranging this yields:

$$R(\vec{k}_{n+1}) = M(\vec{u}_n + \Delta t\vec{k}_{n+1}, \vec{k}_{n+1})+K(\vec{u}_n +\Delta t\vec{k}_{n+1}) - N(t_n + \Delta t)=0$$

Newton iterations must be employed to solve this generally nonlinear system. For Jacobian evaluation, it is useful to write out $M$ as:

$$R(\vec{k}_{n+1}) = \mathbf{M}(\vec{u}_n + \Delta t\vec{k}_{n+1})\vec{k}_{n+1}+K(\vec{u}_n +\Delta t\vec{k}_{n+1}) - N(t_n + \Delta t)=0$$

Then the Jacobian of $R$ is then given by:

$$\dfrac{\partial R(\vec{k}_{n+1})}{\partial \vec{k}}=\mathbf{M}_{n+1} + \Delta t\left.\dfrac{\partial \mathbf{M}}{\partial \vec{u}}\right|_{n+1}\vec{k}_{n+1} + \Delta t \left.\dfrac{\partial K}{\partial \vec{u}}\right|_{n+1} +N(t_n+\Delta t)$$

To be clear, note that, because of the implicit dependence of the $\vec{u}$ argument in $M(\vec{u},\vec{k})$ on $\vec{k}$ here, it must be accounted for in Jacobian calculation, as

$$M(\vec{u}(\vec{k}),\vec{k})=\mathbf{M}(\vec{u}(\vec{k}))\vec{k}$$

$$\dfrac{\partial M}{\partial \vec{k}}=\mathbf{M}(\vec{u}(\vec{k})) + \dfrac{\partial \mathbf{M}}{\partial \vec{u}} \dfrac{\partial \vec{u}}{\partial \vec{k}} $$

and $\partial\vec{u}/\partial\vec{k}=\Delta t$. This means that $\Delta t$ must be given to any nonlinear form integrator that implements it (or must be accounted for if these terms are manually included/accounted for in the $R$ operator).

## Approach #2

### Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. We assume isotropic density $\rho$, thermal conductivity $k$, and specific heat $C$ but allow for them to vary as a function of temperature: $\rho=\rho(u)$, $k=k(u)$, and $C=C(u)$
2. No heat generation

This yields:

$$ \dfrac{\partial u}{\partial t} - \dfrac{1}{\rho C}\nabla \cdot (k \nabla u) = 0$$

We seek a weak formulation of this allowing for nonhomogeneous Neumann and Dirichlet boundary conditions.

Given that one applies the appropriate "lifting" of the solution, we can assume homogeneous Dirichlet BCs in the derivation such that:

$$u=0 \text{ on } \partial \Omega_D$$
$$\hat{n} \cdot (k \nabla u)=g \text{ on } \partial \Omega_N \text{ given } g : \partial \Omega_N \times \mathbb{R}\rightarrow \mathbb{R}$$

where $\Omega$ is our domain and $\partial \Omega$ is the boundary of the domain. Note that $g$ takes in both a location on the boundary and is also allowed to take in a time value that is real.

So: we seek a solution that at a given point in time $u \in H^1_{\partial \Omega_D}(\Omega)$ where $H^1_{\partial \Omega_D}(\Omega) = \{u \in H^1(\Omega) : u=0 \text{ on } \partial \Omega_D \}$. Using conventional notation, the solution space is $\mathcal{V}=H^1_{\partial \Omega_D}(\Omega)$. The test space is chosen to be the same as the solution space.

Multiply the heat equation by test function $v \in H^1_{\partial \Omega_D}(\Omega)$ and integrate:

$$\int_\Omega\dfrac{\partial u}{\partial t}vd\vec{x} - \int_\Omega \dfrac{1}{\rho C}\nabla \cdot (k \nabla u)vd\vec{x} = 0$$

We can then apply the following Green formula to the second term:

$$\int_\Omega \dfrac{1}{\rho C}\nabla \cdot (k \nabla u)vd\vec{x} = -\int_\Omega \dfrac{k}{\rho C}(\nabla v) \cdot (\nabla u) d\vec{x} + \int_{\partial \Omega} \dfrac{1}{\rho C}\hat{n} \cdot (k \nabla u)v d\vec{x}$$

Note from earlier that $ g = \hat{n} \cdot (k \nabla u)$, a specified heat flux on Neumann boundaries. We allow for this heat flux to vary in time as well as in space: $g=g(x,y,z,t)$

Now we may apply a finite element approximation to this weak formulation. This is done assuming $H^1$-conforming nodal elements. An approximate solution $u_h\in\mathcal{V}_h \sub \mathcal{V}$ is sought, where $\mathcal{V}_h$ is a finite-dimensional approximation space, which is a subset of the infinitely-dimensional solution space. Given a basis $\{\phi_j\}$ for the approximation space ($\text{span}\{\phi_j\}=\mathcal{V}_h$) and assuming that the approximate test function $v_h \in \mathcal{V}_h$, we write:
 
$$u(\vec{x},t) \approx u_h(\vec{x},t)= \sum^{N_{\mathcal{V}_h}}_j u_j(t)\phi_j(\vec{x})$$
$$v(\vec{x},t) \approx v_h(\vec{x},t)=\sum^{N_{\mathcal{V}_h}}_j v_j(t)\phi_j(\vec{x})$$


 where the dual basis $u_j$ is chosen to represent $N_{\mathcal{V}_h}$ nodal degrees of freedom in which at a given point in time: $\{u_j\}\in	\mathbb{R}^{N_{\mathcal{V}_h}}$. Going forward for simplification, all summations are presumed to be over ${N_{\mathcal{V}_h}}$ degrees of freedom: $\sum^{N_{\mathcal{V}_h}}_j=\sum_j$.

 Plugging the above in:


$$\int_\Omega \left(\sum_j \dfrac{d u_j(t)}{d t} \phi_j\right) \left(\sum_i v_i(t)\phi_i\right) d\vec{x} + \int_\Omega \left( \dfrac{k(u_h)}{\rho(u_h) C(u_h)}\sum_j u_j(t)\nabla \phi_j\right)\cdot \left(\sum_i v_i(t)\nabla \phi_i\right) d\vec{x} = \int_{\partial \Omega} \dfrac{g}{\rho(u_h) C(u_h)}\left(\sum_i v_i(t) \phi_i \right)d\vec{x}$$


Take out the summations and any non-integrated terms from the integrals:

$$\sum_j \sum_i v_i \left( \int_\Omega \phi_i \phi_j d\vec{x}\right)\dfrac{d u_j}{d t} + \sum_i \sum_j v_i\left( \int_\Omega \dfrac{k(u_h)}{\rho(u_h) C(u_h)}(\nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x}\right)u_j = \sum_i v_i\left( \int_{\partial \Omega} \dfrac{g}{\rho(u_h) C(u_h)}\phi_i d\vec{x}\right)$$

Now we may define:

$$\mathbf{M} = \left[M_{ij}\right]=\left[\int_\Omega \phi_i \phi_j d\vec{x}\right] = \text{Mass Matrix}$$

$$ \mathbf{K}= \left[\kappa_{ij}\right] = \left[\int_\Omega \dfrac{k}{\rho C}(\nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} \right] = \text{Stiffness Matrix}$$

$$\vec{N}=\left[\int_{\partial \Omega} \dfrac{g(\vec{x},t)}{\rho C}\phi_i d\vec{x}\right]=\text{Neumann Term}$$

Thus the equation becomes:

$$\vec{v}^T\mathbf{M}\dfrac{d \vec{u}}{d t} + \vec{v}^T \mathbf{K}\vec{u} = \vec{v}^T \vec{N}$$

$$\vec{v}^T\left(\mathbf{M}\dfrac{d \vec{u}}{d t} + \mathbf{K}\vec{u} - \vec{N}\right)=0$$

For this inner product to be zero $\forall\vec{v}$, the vector in parentheses must be identically zero. So the equation of interest is:
$$
\mathbf{M}\dfrac{d \vec{u}}{d t} + \mathbf{K}\vec{u} - \vec{N}=\vec{0}$$

or

$$M_{ij}\dfrac{du_j}{dt} = -\kappa_{ij}u_j + N_i$$

For linear cases (where all material properties are uniform) the above form with matrices holds. However, generally $\mathbf{K}=\mathbf{K}(\vec{u})$ and $\vec{N}=\vec{N}(\vec{u}, t)$. Thus, to account for these cases, thr equation is written as:

$$\mathbf{M}\dfrac{d\vec{u}}{dt} = A(\vec{u}, t)=-\kappa(\vec{u}) + N(\vec{u},t)$$


where $\kappa$ and $N$, and subsequently $A$, are **nonlinear** operators


### Time-Integration

JOTS makes usage of MFEM's layer-of-abstraction `TimeDependentOperator`.

#### Explicit

For explicit time integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_n$$

where

$$\left.\dfrac{d\vec{u}}{dt}\right|_n=f(\vec{u}_n, t_n)=\mathbf{M}^{-1}A(\vec{u}_n,t_n)$$

#### Implicit

For implicit time-integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_{n+1}$$

where

$$\mathbf{M}\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}=A(\vec{u}_{n+1}, t_{n+1})$$

Plugging in the above equation for $\vec{u}_{n+1}$ and $t_{n+1}=t_n + \Delta t$, and writing $\vec{k}_{n+1}=\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}$

$$\mathbf{M}\vec{k}_{n+1}=A\left(\vec{u}_n + \Delta t \vec{k}_{n+1},t_n + \Delta t\right)=-\kappa(\vec{u}_n + \Delta t\vec{k}_{n+1}) + N(\vec{u}_n + \Delta t \vec{k}_{n+1}, t_n + \Delta t)$$

Rearranging this yields:

$$R(\vec{k}_{n+1}) = \mathbf{M}\vec{k}_{n+1} - A(\vec{u}_n + \Delta t \vec{k}_{n+1},t_n + \Delta t) =0$$

Newton iterations must now be employed to solve this generally nonlinear system. The Jacobian of $R$ is given by:

$$\dfrac{\partial R(\vec{k}_{n+1})}{\partial \vec{k}}=\mathbf{M} - \Delta t \left.\dfrac{\partial A}{\partial \vec{u}}\right|_{n+1}=\mathbf{M}+\Delta t\left.\dfrac{\partial \kappa}{\partial \vec{u}}\right|_{n+1}-\Delta t\left.\dfrac{\partial N}{\partial \vec{u}}\right|_{n+1}$$

## Approximations for Non-Constant Dirichlet + Neumann BCs


Nonhomogeneous time-varying Dirichlet boundary conditions involve fixing some values of $\dfrac{d \vec{u}}{dt}$ **in time**, $\dfrac{d\vec{u}_D}{dt}$, since values of $u$ on $\partial \Omega_D$ may change in time. Currently, JOTS assumes that $\dfrac{d\vec{u}_D}{dt}\approx 0$. This is a source of error as solutions to the derivatives of the free DOFs depend on derivatives of the essential DOFs (See [here](https://github.com/mfem/mfem/issues/1720)).

A similar issue arises for Neumann boundary conditions as JOTS currently assumes that $\vec{N}(t_{n+1})\approx\vec{N}(t_n)$ for implicit time-integration. This is another source of error

Some boundary conditions have an analytical expression, so $\dfrac{d\vec{u}_D}{dt}$ (at $t_n$ if explicit-time, $t_{n+1}$ if implicit-time) and/or $\vec{N}(t_{n+1})$ can be plugged in, but for preCICE boundary conditions, this cannot be done exactly. Approximations may be implemented in the future, such as using a backward differencing for preCICE boundary conditions.

However, for now to maintain consistency, the above assumptions are applied uniformly across JOTS. As $\Delta t \rightarrow0$, the error in the assumptions above also approach zero.

# Steady Heat Transfer

## Weak Formulation Derivation
The following assumptions are made:
1. Isotropic thermal conductivity $k$ allowed to vary as a function of temperature: $k=k(u)$
2. No heat generation

With this, the steady heat equation to be solved is:
$$\nabla \cdot (k \nabla u) = 0$$

As done previously, a weak formulation of this problem is sought. Once again, multiply the heat equation by test function $v$, integrate:

$$\int_\Omega \nabla \cdot (k \nabla u)vd\vec{x} = 0$$

Apply the same Green formula as above to get:

$$\int_\Omega \nabla v \cdot (k\nabla u) d\vec{x} = \int_{\partial \Omega} \hat{n} \cdot (k \nabla u) vd\vec{x}$$

where as before $ g = \hat{n} \cdot (k \nabla u)$, a specified heat flux on Neumann boundaries that *no longer depends on time*.

We no longer presume any dependence on time for the nodal degrees of freedom:

$$u(\vec{x}) \approx u_h(\vec{x})= \sum^{N_{\mathcal{V}_h}}_j u_j\phi_j(\vec{x})$$
$$v(\vec{x}) \approx v_h(\vec{x})=\sum^{N_{\mathcal{V}_h}}_j v_j\phi_j(\vec{x})$$

Plugging in this finite element approximation to the weak formulation gives:

$$\int_\Omega \left( k(u_h)\sum_j u_j\nabla \phi_j\right)\cdot \left(\sum_i v_i\nabla \phi_i\right) d\vec{x} = \int_{\partial \Omega} g\left(\sum_i v_i \phi_i \right)d\vec{x}$$

Rearranging terms:

$$\sum_i \sum_j v_i\left( \int_\Omega (\nabla \phi_i) \cdot (k(u_h) \nabla \phi_j) d\vec{x}\right)u_j = \sum_i v_i\left( \int_{\partial \Omega} g\phi_i d\vec{x}\right)$$

Thus using the same naming conventions as before:

$$\vec{v}^T \mathbf{K}\vec{u} = \vec{v}^T \vec{N}$$

And for this inner product to be identically zero for any test function, we get:

$$\mathbf{K}\vec{u} = \vec{N}$$

where generally this is a nonlinear problem as $\mathbf{K}=\mathbf{K}(\vec{u})$, so:

$$K(\vec{u})=\vec{N}$$

# JOTS Nonlinear Form Integrators

To solve using Newton-Raphson iterations in MFEM using $H^1$-continuous finite elements, new `NonlinearFormIntegrator`'s must be created with member functions `NonlinearFormIntegrator::AssembleElementVector` and `NonlinearFormIntegrator::AssembleElementGrad` implemented. The specific outputs of those functions are shown below for a given LHS nonlinear operator $F=F(\vec{u})$. Note that each operator can be combined into a single LHS operator (as is done for unsteady implicit time-integration). These are used in Newton-Raphson iterations as

$$F(\vec{u}^{k+1}) \approx F(\vec{u}^k) + \dfrac{\partial F(\vec{u}^k)}{\partial \vec{u}}(\vec{u}^{k+1} - \vec{u}^k)$$

$$\rightarrow\vec{u}^{k+1} = \vec{u}^k - \left[\dfrac{\partial F(\vec{u}^k)}{\partial \vec{u}}\right]^{-1}(F(\vec{u}^k) - N_i)$$

For clarity of notation, the gradient $\dfrac{\partial \phi_i}{\partial x_j}$ is now written as $[\nabla \phi]_{ij}$. This was implicitly presumed in the above derivations (that the inner products of matrices are taken), not explicitly included for simplicity. Note that subscript $h$ is not used as an index, but as the finite element approximate solution.

Specifically within JOTS, a derived `NewtonSolver` was created that accepts `MaterialProperty`'s to be updated every Newton iteration within `NewtonSolver::ProcessNewState`. This allows for optimal generality in creation of NLFI's.

## Nonlinear Diffusion Integrator

### `AssembleElementVector`

$K(\vec{u})=\mathbf{K}_{ik}u_k=\displaystyle\int_{\Omega_e} \lambda(u_h)[\nabla\phi]_{il}[\nabla\phi]_{lk} u_kd\vec{x}$

For this, `DiffusionIntegrator` is simply used with the coefficient $\lambda(u_h)$ set as its `Coefficient`. The action of the operator on $u_k$ is then computed with `DiffusionIntegrator::AssembleElementVector`.


### `AssembleElementGrad`

$\dfrac{\partial K(\vec{u})}{\partial \vec{u}}= \dfrac{\partial}{\partial u_j}\left(\mathbf{K}_{ik}u_k\right) = \dfrac{\partial \mathbf{K}_{ik}}{\partial u_j}u_k + \mathbf{K}_{ik}\delta_{kj}= \displaystyle\int_{\Omega_e} \lambda'(u_h)[\nabla\phi]_{il}[\nabla\phi]_{lk}u_k \phi_j d\vec{x} + \displaystyle\int_{\Omega_e} \lambda(u_h)[\nabla\phi]_{il} [\nabla\phi]_{lj} d\vec{x}$

For the first term, a `MixedScalarWeakDivergenceIntegrator::AssembleElementMatrix` with a `ScalarVectorProductCoefficient` is used, multiplying the $\lambda'(u_h)$ coefficient with the gradient of the solution represented as a `GradientGridFunctionCoefficient`. For the second term, the same `DiffusionIntegrator` from before is used and `DiffusionIntegrator::AssembleElementMatrix` is called.


## Nonlinear Neumann, $N$

### `AssembleElementVector`
$N(\vec{u}, t)=\displaystyle \int_{\partial \Omega_e} \lambda(\vec{u},\vec{x},t)\phi_i d\vec{x}$


### `AssembleElementGrad`
$\dfrac{\partial N}{\partial u_j}=\displaystyle \int_{\partial \Omega_e}\dfrac{\partial\lambda(\vec{u},\vec{x},t)}{\partial \vec{u}}\phi_i \phi_j d\vec{x}$

# Notes:

- Presently:
    - Times that have a difference less than $10^{-14}$ are presumed equal
    - Restart files output data to 15 decimal points