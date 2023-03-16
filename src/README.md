# Dev Notes


## Understanding the discretization

It is quite complex. See this website for the derivation: https://en.wikiversity.org/wiki/Finite_elements/Solution_of_heat_equation

Or you can look at this much simpler derivation by the one and only ChatGPT:

Starting from the continuous form of the unsteady thermal conduction equation:

$$\rho c_p \frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = f$$

where $u$ is the temperature, $\rho$ is the density, $c_p$ is the specific heat, $k$ is the thermal conductivity, and $f$ is a volumetric heat source/sink term.

We begin by approximating the temperature $u$ using a finite element basis expansion:

$$u(x) \approx \sum_i u_i \phi_i(x)$$

where $u_i$ are the coefficients of the expansion and $\phi_i(x)$ are the finite element basis functions.

Substituting this expression for $u$ into the continuous equation and multiplying by a test function $\phi_j(x)$, we obtain:

$$\int_\Omega \left(\rho c_p \frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u)\right) \phi_j(x) \ dx = \int_\Omega f \phi_j(x) \ dx$$

Using integration by parts to eliminate the second-order derivative term, we get:

$$\int_\Omega \rho c_p \frac{\partial u}{\partial t} \phi_j(x) \ dx + \int_\Omega k \nabla u \cdot \nabla \phi_j(x) \ dx - \int_{\partial \Omega} k \frac{\partial u}{\partial n} \phi_j(x) \ ds = \int_\Omega f \phi_j(x) \ dx$$

where $\partial \Omega$ denotes the boundary of the domain $\Omega$ and $\frac{\partial u}{\partial n}$ is the normal derivative of $u$ on the boundary.

Assuming that the boundary conditions satisfy the homogeneous Neumann condition ($\frac{\partial u}{\partial n} = 0$), the boundary term vanishes and we are left with:

$$\int_\Omega \rho c_p \frac{\partial u}{\partial t} \phi_j(x) \ dx + \int_\Omega k \nabla u \cdot \nabla \phi_j(x) \ dx = \int_\Omega f \phi_j(x) \ dx$$

Using the Galerkin finite element method, we approximate the temperature $u$ and test function $\phi_j(x)$ using the same basis functions and restrict the integrals to the finite element subdomains:

$$\sum_i \left(\rho c_p \frac{\partial u_i}{\partial t} + \int_\Omega k \nabla \phi_i(x) \cdot \nabla \phi_j(x) \ dx \ u_i\right) = \int_\Omega f \phi_j(x) \ dx$$

This can be written more compactly as:

$$\sum_i (M_{ij} \frac{du_i}{dt} + K_{ij} u_i) = F_j$$

where $M_{ij}$ is the mass matrix and $K_{ij}$ is the stiffness matrix, given by:

$$M_{ij} = \int_\Omega \rho c_p \phi_i(x) \phi_j(x) \ dx, \ \ K_{ij} = \int_\Omega k \nabla \phi_i(x) \cdot \nabla \phi_j(x) \ dx$$

The final step is to apply the inverse of the mass matrix $M$ to both sides of the equation, giving:

$$M^{-1} \sum_i (M_{ij} \frac{du_i}{dt} + K_{ij} u_i) = M^{-1} F_j$$

Using the distributive property of matrix multiplication, we can simplify this equation to:

$$\frac{d}{dt} \sum_i (M^{-1} M_{ij} u_i) + \sum_i (M^{-1} K_{ij} u_i) = M^{-1} F_j$$

Recognizing that the term $\sum_i (M^{-1} M_{ij} u_i)$ is just the finite element approximation of the temperature $u$ at the point $x_j$, we can rewrite this equation as:

$$\frac{d}{dt} u_j + \sum_i (M^{-1} K_{ij} u_i) = (M^{-1} F)_j$$

This is the desired equation for the finite element method after spatial discretization:

$$\frac{d}{dt} u = M^{-1}(-Ku) + M^{-1} F$$

where $u$ is the temperature vector, $M$ is the mass matrix, $K$ is the stiffness matrix, and $F$ is the load vector. This equation can be solved numerically using time integration methods such as the backward Euler method or the Crank-Nicolson method.

Note that for our case, we are assuming $F=0$