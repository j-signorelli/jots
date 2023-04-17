# Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. Density $\rho$ and specific heat $C_p$ are uniform (constant in space and time)
2. We assume isotropic thermal conductivity $k$ but allow for it to vary as a function of temperature: $k=k(T(x,y,z,t))$
3. No heat generation

This yields:

$$\rho c_p \dfrac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = 0$$

We seek a weak formulation of this allowing for nonhomogeneous Neumann and Dirichlet boundary conditions.

Given that one applies the appropriate "lifting" of the solution, we can assume homogeneous Dirichlet BCs in the derivation such that:

$$T=0 \text{ on } \partial \Omega_D$$
$$\hat{n} \cdot (k \nabla T)=g \text{ on } \partial \Omega_N \text{ given } g : \partial \Omega_N \rightarrow \mathbb{R}$$

where $\Omega$ is our domain and $\partial \Omega$ is the boundary of the domain.

So: we seek a solution $T \isin H^1_{\partial \Omega_D}(\Omega)$ where $H^1_{\partial \Omega_D}(\Omega) = \{T \isin H^1(\Omega) : T=0 \text{ on } \partial \Omega_D \}$.

Multiply the heat equation by test function $v \isin H^1_{\partial \Omega_D}(\Omega)$ and integrate:

$$\int_\Omega \rho c_p \dfrac{\partial T}{\partial t}vd\vec{x} - \int_\Omega \nabla \cdot (k \nabla T)vd\vec{x} = 0$$

We can then apply the following Green formula to the second term:

$$\int_\Omega \nabla \cdot (k \nabla T)vd\vec{x} = -\int_\Omega \nabla v \cdot (k\nabla T) d\vec{x} + \int_{\partial \Omega} \hat{n} \cdot (k \nabla T) d\vec{x}$$

Note from earlier that $ g = \hat{n} \cdot (k \nabla T)$, a specified heat flux on Neumann boundaries.

Now we may apply the finite element approximation to this weak formulation, assuming $H^1$-conforming nodal elements:

$$T(x,y,z,t) = \sum_j T_j(t)\phi_j(x,y,z)$$
$$v(x,y,z) = \sum_j v_j(t)\phi_j(x,y,z)$$


where $\phi_j(x)$ are the finite element basis functions and $T_j$ are nodal DOFs. Plugging this in:

$$\sum_i \sum_j \left( \int_\Omega \rho C_p \dfrac{d T_i}{d t} \phi_i \phi_j d\vec{x} v_j\right) + \sum_i \sum_j \left( \int_\Omega T_i (k \nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} v_j\right) = \sum_j \left( \int_{\partial \Omega} g\phi_j d\vec{x}v_j\right)$$

Note now that we may define:

$$M_{ij}=\int_\Omega \rho C_p \phi_i \phi_j d\vec{x} = \text{Mass Matrix}$$

$$K_{ij} = \int_\Omega (k \nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} = \text{Stiffness Matrix}$$

and then rewrite:

$$M_{ij}\dfrac{dT_j}{dt} = -K_{ij}T_j + \int_{\partial \Omega} g\phi_i d\vec{x}$$

Finally, we can rewrite this as:

$$\dfrac{dT_i}{dt} = M_{ik}^{-1} \left( -K_{kj}T_j + \int_{\partial \Omega} g\phi_k d\vec{x}\right)$$

Note that the boundary integral term is a linear form that must be added to enforce Neumann BCs.

# Nonhomogeneous Dirichlet BCs

To enforce nonhomogeneous Dirichlet BCs, given temperatures at essential DOFs $T_e$ are enforced, with $\dfrac{\partial T_e}{\partial t} = 0$ as of now.

It may be straightforward to update this for any Dirichlet BCs that change in time, where one may apply a backward differencing if desired, but calculations at start given an initial condition must be carefully considered.

# Time-Integration

For explicit time-integration schemes, it is simply assumed that the $T_j$ on the RHS multiplied by $K$ is the temperature at the previous timestep. For implicit time-integration schemes, $T_j$ on the RHS must be the temperature at the next timestep, so the equation becomes:

$$\dfrac{dT_i}{dt} = M_{ik}^{-1} \left( -K_{kj}(T_j + dt\dfrac{dT_k}{dt}) + \int_{\partial \Omega} g\phi_k d\vec{x}\right)$$

as taking the Taylor Series of $T_i^n$ about the $t_{n+1}$ timestep yields $T_i^n = T_i^{n+1} - dt\dfrac{\partial T^{n+1}_i}{\partial t} + O(dt^2)$ and neglecting HOT yields:  


$$T_i^{n+1} = T_i^n + dt \dfrac{dT^{n+1}_i}{dt}$$


Given this, the new equation to solve becomes:

$$Kdt\dfrac{dT}{dt} + M\dfrac{dT}{dt} = \left( -KT + \int_{\partial \Omega} g\phi d\vec{x}\right)$$

$$\left( Kdt + M\right)\dfrac{dT}{dt} = \left( -KT + \int_{\partial \Omega} g\phi d\vec{x}\right)$$

Note **importantly**: stiffness matrix $K$ is still calculated using the temperature from the previous timestep so, if it changes in time, this is a source of error.