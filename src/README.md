
# Weak Formulation Derivation

Start with the unsteady thermal conduction equation, applying the following assumptions:

1. Density $\rho$ is constant in space and time
2. We assume isotropic thermal conductivity $k$ and specific heat $C$ but allow for them to vary as a function of temperature: $k=k(T)$ and $C=C(T)$
3. No heat generation

This yields:

$$\rho C \dfrac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = 0$$

We seek a weak formulation of this allowing for nonhomogeneous Neumann and Dirichlet boundary conditions.

Given that one applies the appropriate "lifting" of the solution, we can assume homogeneous Dirichlet BCs in the derivation such that:

$$T=0 \text{ on } \partial \Omega_D$$
$$\hat{n} \cdot (k \nabla T)=g \text{ on } \partial \Omega_N \text{ given } g : \partial \Omega_N \times \mathbb{R}\rightarrow \mathbb{R}$$

where $\Omega$ is our domain and $\partial \Omega$ is the boundary of the domain. Note that $g$ takes in both a location on the boundary and is also allowed to take in a time value that is real.

So: we seek a solution $T \in H^1_{\partial \Omega_D}(\Omega)$ where $H^1_{\partial \Omega_D}(\Omega) = \{T \in H^1(\Omega) : T=0 \text{ on } \partial \Omega_D \}$. Using conventional notation, the solution space is $\mathcal{V}=H^1_{\partial \Omega_D}(\Omega)$. The test space is chosen to be the same as the solution space.

Multiply the heat equation by test function $v \in H^1_{\partial \Omega_D}(\Omega)$ and integrate:

$$\int_\Omega \rho C \dfrac{\partial T}{\partial t}vd\vec{x} - \int_\Omega \nabla \cdot (k \nabla T)vd\vec{x} = 0$$

We can then apply the following Green formula to the second term:

$$\int_\Omega \nabla \cdot (k \nabla T)vd\vec{x} = -\int_\Omega \nabla v \cdot (k\nabla T) d\vec{x} + \int_{\partial \Omega} \hat{n} \cdot (k \nabla T) d\vec{x}$$

Note from earlier that $ g = \hat{n} \cdot (k \nabla T)$, a specified heat flux on Neumann boundaries. We allow for this heat flux to vary in time as well as in space: $g=g(x,y,z,t)$

Now we may apply a finite element approximation to this weak formulation. This is done assuming $H^1$-conforming nodal elements. An approximate solution $T_h\in\mathcal{V}_h \sub \mathcal{V}$ is sought, where $\mathcal{V}_h$ is a finite-dimensional approximation space, which is a subset of the infinitely-dimensional solution space. Given a basis $\{\phi_j\}$ for the approximation space ($\text{span}\{\phi_j\}=\mathcal{V}_h$) and assuming that the approximate test function $v_h \in \mathcal{V}_h$, we write:
 
$$T(\vec{x},t) \approx T_h(\vec{x},t)= \sum^{N_{\mathcal{V}_h}}_j T_j(t)\phi_j(\vec{x})$$
$$v(\vec{x},t) \approx v_h(\vec{x},t)=\sum^{N_{\mathcal{V}_h}}_j v_j(t)\phi_j(\vec{x})$$


 where $T_j$ represent $N_{\mathcal{V}_h}$ nodal degrees of freedom in which at a given point in time: $\{T_j\}\in	\mathbb{R}^{N_{\mathcal{V}_h}}$. Going forward for simplification, all summations are presumed to be over ${N_{\mathcal{V}_h}}$ degrees of freedom: $\sum^{N_{\mathcal{V}_h}}_j=\sum_j$. Plugging the above in:


$$\int_\Omega \rho C \left(\sum_j \dfrac{d T_j(t)}{d t} \phi_j\right) \left(\sum_i v_i(t)\phi_i\right) d\vec{x} + \int_\Omega \left( k\sum_j T_j(t)\nabla \phi_j\right)\cdot \left(\sum_i v_i(t)\nabla \phi_i\right) d\vec{x} = \int_{\partial \Omega} g\left(\sum_i v_i(t) \phi_i \right)d\vec{x}$$

Take out the summations and any non-integrated terms from the integrals:

$$\sum_j \sum_i v_i \left( \int_\Omega \rho C \phi_i \phi_j d\vec{x}\right)\dfrac{d T_j}{d t} + \sum_i \sum_j v_i\left( \int_\Omega (\nabla \phi_i) \cdot (k \nabla \phi_j) d\vec{x}\right)T_j = \sum_i v_i\left( \int_{\partial \Omega} g\phi_i d\vec{x}\right)$$

Now we may define:

$$\bold{M} = \left[M_{ij}\right]=\left[\int_\Omega \rho C \phi_i \phi_j d\vec{x}\right] = \text{Mass Matrix}$$

$$\bold{K} = \left[K_{ij}\right] = \left[\int_\Omega (k \nabla \phi_i) \cdot (\nabla \phi_j) d\vec{x} \right] = \text{Stiffness Matrix}$$

$$\vec{N}=\left[\int_{\partial \Omega} g(\vec{x},t)\phi_i d\vec{x}\right]=\text{Neumann Linear Form}$$

^^HOW DO I WRITE DEPENDENCIES ON NODAL TEMPERATURE IN THE ABOVE BILINEAR FORMS?!??

Thus the equation becomes:

$$\vec{v}^T\bold{M}\dfrac{d \vec{T}}{d t} = \vec{v}^T \bold{K}\vec{T} + \vec{v}^T \vec{N}$$

and then rewrite while dividing out :

$$M_{ij}\dfrac{dT_j}{dt} = -K_{ij}T_j + \int_{\partial \Omega} g\phi_i d\vec{x}$$

Finally, we can rewrite this as:

$$\dfrac{dT_i}{dt} = M_{ik}^{-1} \left( -K_{kj}T_j + \int_{\partial \Omega} g\phi_k d\vec{x}\right)$$

Note that the boundary integral term is a linear form that must be added to enforce Neumann BCs.

# Local Linearization

To aid in time-integration, a locally-linearized form of the above equation is sought.

Write the $F(T,t)=M^{-1}(T)$ and Taylor expand $F(T,t)$ and $T(t)$ about a reference point in time $t_n$:

$$F(T,t) = F(T_n,t_n) + \left. \dfrac{\partial F}{\partial T}\right|_n(T-T_n) + \left. \dfrac{\partial F}{\partial t}\right|_n(t-t_n) + \dfrac{1}{2}\left.\dfrac{\partial^2 F}{\partial T^2}\right|_n(T-T_n)^2 + \dfrac{1}{2}\left.\dfrac{\partial^2 F}{\partial t^2}\right|_n(t-t_n)^2 \dfrac{1}{2}\left.\dfrac{\partial^2 F}{\partial T \partial t}\right|_n(T-T_n)(t-t_n)+...$$

and

$$T(t) = T_n + (t-t_n)\left.\dfrac{\partial T}{\partial t}\right|_n + \dfrac{1}{2}\left.\dfrac{\partial^2 T}{\partial t^2}\right|_n+...$$

Note that if $t$ is within $\Delta t$ of $t_n$, then both $(t-t_n)^k$ and $(T-T_n)^k$ are $\mathcal{O}(\Delta t^k)$, so we can write:

$$\dfrac{dT}{dt}=F(T,t)=F(T_n,t_n) + \left. \dfrac{\partial F}{\partial T}\right|_n(T-T_n) + \left. \dfrac{\partial F}{\partial t}\right|_n(t-t_n) + \mathcal{O}(\Delta t^2)$$

Note that the mass matrix $M$ and stiffness matrix $K$ do not explicitly depend on $t$, but it is possible that the Neumann linear form does.

# Unsteady Simulations -- Time-Integration

## Explicit
For explicit time-integration we seek $T^{n+1}=G(T^{n})$. In terms of implicit and explicit parts $F$ and $G$ respectively we have:

$$F(T^{n+1})=\dfrac{dT}{dt}=\dfrac{T^{n+1} - T^n}{\Delta t}$$

and
$$G(T^n)=M^{-1}(T^n) \left( -K(T^n)T^n + \int_{\partial \Omega} g\phi^n d\vec{x}\right)$$

Thus:

$$\dfrac{dT}{dt}^{n+1}=M^{-1}(T^n) \left( -K(T^n)T^n + \int_{\partial \Omega} g\phi^n d\vec{x}\right)$$

## Implicit
Implicit time-integration schemes require solving $H(T^{n+1}, T^n)=0$

Either an implicit-explicit (IMEX) approach or a nonlinear solver must be used.

### IMEX

In the IMEX approach, we seek $T^{n+1}=F(T^{n+1}) + G(T^n)$ where, as before, $F$ and $G$ are the implicit and explicit parts respectively. The implicit part is chosen to include linear terms while the explicit part is chosen to include nonlinear terms.



Linearizing $T_j^{n+1}$ yields the equation:

$$M_{ik}\dfrac{dT_i}{dt} = \left( -K_{kj}(T_j + dt\dfrac{dT_k}{dt}) + \int_{\partial \Omega} g\phi_k d\vec{x}\right)$$

as taking the Taylor Series of $T_i^n$ about the $t_{n+1}$ timestep yields $T_i^n = T_i^{n+1} - dt\dfrac{\partial T^{n+1}_i}{\partial t} + O(dt^2)$ and neglecting HOT yields:  


$$T_i^{n+1} = T_i^n + dt \dfrac{dT^{n+1}_i}{dt}$$


Given this, the new equation to solve becomes:

$$Kdt\dfrac{dT}{dt} + M\dfrac{dT}{dt} = \left( -KT + \int_{\partial \Omega} g\phi d\vec{x}\right)$$

$$\left( Kdt + M\right)\dfrac{dT}{dt} = \left( -KT + \int_{\partial \Omega} g\phi d\vec{x}\right)$$

Note **importantly**: stiffness matrix $K$ is still calculated using the temperature from the previous timestep so, if it changes in time, this is a source of error. This is beneficial as only one matrix must be inverted each timestep.

Alby Questions:
- What $g$ do I use??? $g^{n+1}$ or $g^n$??
- What exactly is being done in the "Linearized approximation"? I want to specifically see all the math written out and understand the names of these things (IMEX?)

### Picard Iterations
Really simple nonlinear solver with awful convergence but easy.

### Newton-Raphson
If we have a Jacobian, then perfect. I think I must provide that. Otherwise we need to calculate a Jacobian which may be expensive.

Then, need to invert the Jacobian.

## Steady-State Simulations
The steady state equation of interest to solve becomes:

$$\nabla \cdot (k \nabla T) = 0$$

If $k=k(T)$, then a nonlinear solver must be used. Either Picard iterations or Newton-Raphson may be used.

# Nonhomogeneous Dirichlet BCs

To enforce nonhomogeneous Dirichlet BCs, given temperatures at essential DOFs $T_e$ are enforced, with $\dfrac{\partial T_e}{\partial t} = 0$ as of now.

It may be straightforward to update this for any Dirichlet BCs that change in time, where one may apply a backward differencing if desired, but calculations at start given an initial condition must be carefully considered.


# Important Assumptions / Limitations

1. For implicit time-integration, $K$ is determined using temperatures at the previous timestep ($K$ is linearized from previous timestep).
2. For essential boundary conditions that vary in time, $\dfrac{dT_i}{dt}$ is assumed to be 0. A backward differencing may be implemented for higher-order accuracy in time but would require restart files that save $t_{n-1}$ data in addition to $t_n$ data.

# Notes:

- Presently:
    - Times that have a difference less than $10^{-14}$ are presumed equal
    - Restart files output data to 15 decimal points