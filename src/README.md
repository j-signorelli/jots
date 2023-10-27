# Unsteady Heat Transfer

## Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. Density $\rho$ is constant in space and time
2. We assume isotropic thermal conductivity $k$ and specific heat $C$ but allow for them to vary as a function of temperature: $k=k(u)$ and $C=C(u)$
3. No heat generation

This yields:

$$\rho C \dfrac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = 0$$

We seek a weak formulation of this allowing for nonhomogeneous Neumann and Dirichlet boundary conditions.

Given that one applies the appropriate "lifting" of the solution, we can assume homogeneous Dirichlet BCs in the derivation such that:

$$u=0 \text{ on } \partial \Omega_D$$
$$\hat{n} \cdot (k \nabla u)=g \text{ on } \partial \Omega_N \text{ given } g : \partial \Omega_N \times \mathbb{R}\rightarrow \mathbb{R}$$

where $\Omega$ is our domain and $\partial \Omega$ is the boundary of the domain. Note that $g$ takes in both a location on the boundary and is also allowed to take in a time value that is real.

So: we seek a solution that at a given point in time $u \in H^1_{\partial \Omega_D}(\Omega)$ where $H^1_{\partial \Omega_D}(\Omega) = \{u \in H^1(\Omega) : u=0 \text{ on } \partial \Omega_D \}$. Using conventional notation, the solution space is $\mathcal{V}=H^1_{\partial \Omega_D}(\Omega)$. The test space is chosen to be the same as the solution space.

Multiply the heat equation by test function $v \in H^1_{\partial \Omega_D}(\Omega)$ and integrate:

$$\int_\Omega \rho C \dfrac{\partial u}{\partial t}vd\vec{x} - \int_\Omega \nabla \cdot (k \nabla u)vd\vec{x} = 0$$

We can then apply the following Green formula to the second term:

$$\int_\Omega \nabla \cdot (k \nabla u)vd\vec{x} = -\int_\Omega \nabla v \cdot (k\nabla u) d\vec{x} + \int_{\partial \Omega} \hat{n} \cdot (k \nabla u)v d\vec{x}$$

Note from earlier that $ g = \hat{n} \cdot (k \nabla u)$, a specified heat flux on Neumann boundaries. We allow for this heat flux to vary in time as well as in space: $g=g(x,y,z,t)$

Now we may apply a finite element approximation to this weak formulation. This is done assuming $H^1$-conforming nodal elements. An approximate solution $u_h\in\mathcal{V}_h \sub \mathcal{V}$ is sought, where $\mathcal{V}_h$ is a finite-dimensional approximation space, which is a subset of the infinitely-dimensional solution space. Given a basis $\{\phi_j\}$ for the approximation space ($\text{span}\{\phi_j\}=\mathcal{V}_h$) and assuming that the approximate test function $v_h \in \mathcal{V}_h$, we write:
 
$$u(\vec{x},t) \approx u_h(\vec{x},t)= \sum^{N_{\mathcal{V}_h}}_j u_j(t)\phi_j(\vec{x})$$
$$v(\vec{x},t) \approx v_h(\vec{x},t)=\sum^{N_{\mathcal{V}_h}}_j v_j(t)\phi_j(\vec{x})$$


 where the dual basis $u_j$ is chosen to represent $N_{\mathcal{V}_h}$ nodal degrees of freedom in which at a given point in time: $\{u_j\}\in	\mathbb{R}^{N_{\mathcal{V}_h}}$. Going forward for simplification, all summations are presumed to be over ${N_{\mathcal{V}_h}}$ degrees of freedom: $\sum^{N_{\mathcal{V}_h}}_j=\sum_j$.

 Plugging the above in:


$$\int_\Omega \rho C(u_h) \left(\sum_j \dfrac{d u_j(t)}{d t} \phi_j\right) \left(\sum_i v_i(t)\phi_i\right) d\vec{x} + \int_\Omega \left( k(u_h)\sum_j u_j(t)\nabla \phi_j\right)\cdot \left(\sum_i v_i(t)\nabla \phi_i\right) d\vec{x} = \int_{\partial \Omega} g\left(\sum_i v_i(t) \phi_i \right)d\vec{x}$$


Take out the summations and any non-integrated terms from the integrals:

$$\sum_j \sum_i v_i \left( \int_\Omega \rho C(u_h) \phi_i \phi_j d\vec{x}\right)\dfrac{d u_j}{d t} + \sum_i \sum_j v_i\left( \int_\Omega (\nabla \phi_i) \cdot (k(u_h) \nabla \phi_j) d\vec{x}\right)u_j = \sum_i v_i\left( \int_{\partial \Omega} g\phi_i d\vec{x}\right)$$

Now we may define:

$$\bold{M} = \left[M_{ij}\right]=\left[\int_\Omega \rho C \phi_i \phi_j d\vec{x}\right] = \text{Mass Matrix}$$

$$\bold{K} = \left[K_{ij}\right] = \left[\int_\Omega (k \nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} \right] = \text{Stiffness Matrix}$$

$$\vec{N}=\left[\int_{\partial \Omega} g(\vec{x},t)\phi_i d\vec{x}\right]=\text{Neumann Linear Form}$$

Thus the equation becomes:

$$\vec{v}^T\bold{M}\dfrac{d \vec{u}}{d t} + \vec{v}^T \bold{K}\vec{u} = \vec{v}^T \vec{N}$$

$$\vec{v}^T\left(\bold{M}\dfrac{d \vec{u}}{d t} + \bold{K}\vec{u} - \vec{N}\right)=0$$

For this inner product to be zero $\forall\vec{v}$, the vector in parentheses must be identically zero. So the equation of interest is:
$$
\bold{M}\dfrac{d \vec{u}}{d t} + \bold{K}\vec{u} - \vec{N}=\vec{0}$$

or

$$M_{ij}\dfrac{du_j}{dt} = -K_{ij}u_j + N_i$$

Finally, we can rewrite this as:

$$\dfrac{du_i}{dt} = M_{ik}^{-1} \left( -K_{kj}u_j + N_k\right)=-M_{ik}^{-1}K_{kj}u_j + M_{ik}^{-1}N_k$$

Note that generally $\bold{M}=\bold{M}(\vec{u})$,  $\bold{K}=\bold{K}(\vec{u})$, and $\vec{N}=\vec{N}(t)$

## Time-Integration

JOTS makes usage of MFEM's layer-of-abstraction `TimeDependentOperator`. Note that $k=\dfrac{d\vec{u}}{dt}$

### Explicit

For explicit time integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_n$$

where

$$\left.\dfrac{d\vec{u}}{dt}\right|_n=f(\vec{u}_n, t_n)=\bold{M}^{-1}(\vec{u}_n)\left[\bold{K}(\vec{u}_n)\vec{u}_{n} + \vec{N}(t_n)\right]$$

### Implicit

For implicit time-integration methods,

$$\vec{u}_{n+1} = \vec{u}_n + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_{n+1}$$

where

$$\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}=f(\vec{u}_{n+1}, t_{n+1})=\bold{M}^{-1}(\vec{u}_{n+1})\left[\bold{K}(\vec{u}_{n+1})\vec{u}_{n+1} + \vec{N}(t_{n+1})\right]$$

Plugging in the above equation for $\vec{u}_{n+1}$:

$$\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}=f(\vec{u}_{n+1}, t_{n+1})=\bold{M}^{-1}(\vec{u}_{n+1})\left[\bold{K}(\vec{u}_{n+1})\left(\vec{u}_{n} + \Delta t\left.\dfrac{d \vec{u}}{dt}\right|_{n+1}\right) + \vec{N}(t_{n+1})\right]$$

Rearranging this yields:

$$\left[\bold{M}(\vec{u}_{n+1})+\Delta t \bold{K}(\vec{u}_{n+1})\right]\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}=\bold{K}(\vec{u}_{n+1})\vec{u}_{n}+\vec{N}(t_{n+1})$$

Thus:

$$\left.\dfrac{d\vec{u}}{dt}\right|_{n+1}=\left[\bold{M}(\vec{u}_{n+1})+\Delta t \bold{K}(\vec{u}_{n+1})\right]^{-1}\left[\bold{K}(\vec{u}_{n+1})\vec{u}_{n}+\vec{N}(t_{n+1})\right]$$

Given the nonlinearities, JOTS presently only has one way to deal with them by considering the mass and stiffness matrices' Taylor expansion taken about the previous timestep:

$$\bold{M}^{-1}(\vec{u}_{n+1})=\bold{M}^{-1}(\vec{u}_{n}) + \left.\dfrac{\partial \bold{M}^{-1}}{\partial \vec{u}}\right|_n\left[ \vec{u}_{n+1} - \vec{u}_{n}\right] + \dfrac{1}{2}\left.\dfrac{\partial^2 \bold{M}^{-1}}{\partial \vec{u}^2}\right|_n\left[ \vec{u}_{n+1} - \vec{u}_{n}\right]^2 + \dotsb$$

$$\bold{K}(\vec{u}_{n+1})=\bold{K}(\vec{u}_{n}) + \left.\dfrac{\partial \bold{K}}{\partial \vec{u}}\right|_n\left[ \vec{u}_{n+1} - \vec{u}_{n}\right] + \dfrac{1}{2}\left.\dfrac{\partial^2 \bold{K}}{\partial \vec{u}^2}\right|_n\left[ \vec{u}_{n+1} - \vec{u}_{n}\right]^2 + \dotsb$$

It is currently presumed that:

$$\bold{M}^{-1}(\vec{u}_{n+1})=\bold{M}^{-1}(\vec{u}_{n})$$

and

$$\bold{K}(\vec{u}_{n+1})=\bold{K}(\vec{u}_{n}).$$


There are plans to implement nonlinear iterative solvers for these terms and this section will be expanded accordingly.

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

$$\vec{v}^T \bold{K}\vec{u} = \vec{v}^T \vec{N}$$

And for this inner product to be identically zero for any test function, we get:

$$\bold{K}\vec{u} = \vec{N}$$

where generally this is a nonlinear problem as $\bold{K}=\bold{K}(\vec{u})$.

# Notes:

- Presently:
    - Times that have a difference less than $10^{-14}$ are presumed equal
    - Restart files output data to 15 decimal points