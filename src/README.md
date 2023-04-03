# Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. Density $\rho$ and specific heat $C_p$ are uniform (constant in space and time)
2. We allow for thermal conductivity $k$ to vary as a function of temperature: $k=k(T(x,y,z,t))$
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

$$M_{ij}\dfrac{dT_j}{dt} + K_{ij}T_j = \int_{\partial \Omega} g\phi_j d\vec{x}$$

Finally, we can rewrite this as:

$$\dfrac{dT_j}{dt} = M_{ij}^{-1} \left( -K_{ij}T_j + \int_{\partial \Omega} g\phi_j d\vec{x}\right)$$

Note that the boundary integral term is a linear form that must be added to enforce Neumann BCs.