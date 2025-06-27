# Chapter 4: The Finite Volume library

[Back to the table of contents](./0_start.md)
or [Previous chapter: Making a mesh](./3_mesh.md)

## Finite volume discretization in Briscola

As was mentioned in [Chapter 3](./3_mesh.md), Briscola uses the Finite Volume
Method (FVM) for spatial discretization. In the FVM, the domain is discretized
into multiple control volumes (i.e., mesh cells), and governing equations are
integrated over each control volume to obtain their integral form. Fluxes of
conserved quantities are defined at cell faces, such that they can be integrated
over the control volumes by summation over all of the faces of each control
volume (using Gauss's theorem). This also ensures consistency of fluxes between
neighboring control volumes. In Briscola, the classes and functions related to
the FVM are grouped into the `briscolaFiniteVolume` library.

## Boundary conditions

In Briscola, boundary conditions can broadly be grouped into two categories:
standard boundary conditions and immersed boundary conditions. Standard boundary
conditions are defined at domain boundary patches and therefore always coincide
exactly with cell faces. Immersed boundary conditions on the other hand are
internal to the flow domain, and can be set on immersed boundaries which may
arbitrarily cross cells and cell faces. The immersed boundary method and
immersed boundary conditions are discussed in more detail in [Chapter
6](./6_ibm.md).

For standard boundary conditions, Briscola uses the ghost-cell method [Harlow,
F. H., & Welch, J. E. (1965) Phys. Fluids 8(12) : 2182]. In the ghost-cell
method, each processor's partial mesh is extended with one layer of padding
cells on all of its patches. These so-called ghost-cells are not part of the
flow domain, and the values of flow quantities are not solved for in these
cells. Instead, their values are set such that, according to the finite volume
discretization across the domain boundary, the desired boundary condition is
enforced at the domain boundary. The main advantage of the ghost-cell method is
that the same discretization can be used throughout the entire mesh, including
for boundary-adjacent cells. As Briscola uses a compact discretization stencil,
a single layer of ghost-cells is required. However, for higher order
discretizations with larger stencils, the number of ghost-cell layers could
simply be increased, although this is not yet implemented. The following
boundary condition types are available in Briscola:

* Dirichlet boundary condition: a fixed boundary value $\phi_b$ needs to be
specified. The ghost cell value $\phi_g$ is set as $\phi_g = 2\phi_b - \phi_f$,
with $\phi_f$ the value in the boundary-adjacent cell inside the flow domain.
* Neumann boundary condition: a fixed boundary-normal gradient $g_b$ needs to be
specified. The ghost cell value $\phi_g$ is set as $\phi_g = \Delta x \cdot
g_b+\phi_f$, with $\Delta x$ the cell spacing, and $\phi_f$ the value in the
boundary-adjacent cell inside the flow domain.
* Empty boundary condition: For two-dimensional meshes, the empty boundary
condition must be used in the unused third spatial direction. This boundary
condition is a specialization of the Neumann boundary condition with a
zero-gradient.
* No-slip boundary condition: This boundary condition is a specialization of the
Dirichlet boundary condition with a boundary value of zero.
* Slip boundary condition: This boundary condition is meant for the velocity at
a slip wall. It sets the boundary-normal velocity and the wall shear to zero.
* Outflow boundary condition: This boundary condition is also a specialization
of the Neumann boundary condition with a zero-gradient. However, for staggered
quantities, it also checks for and prevents inflow. This is needed for numerical
stability of staggered flow solvers in certain cases.
* Parallel boundary condition: This boundary condition is generated
automatically when a mesh is decomposed across multiple processors. It sets the
ghost-cell values to be equal to the neighboring processor's boundary-adjacent
cell values to ensure the continuity of the solution across multiple processors.
Communication between neighboring processors is handled with OpenFOAM's
`PStream` library, which itself uses MPI.
* Periodic boundary condition: This boundary condition is a specialization of
the parallel boundary condition, but where communication occurs not between
neighboring processors but between processors adjacent to a pair of periodic
patches. This can also be one and the same processor if there is no parallel
decomposition in the periodic direction.
* Dummy boundary condition: This boundary condition is used as a dummy for
fields which don't need any particular boundary conditions. Correcting this
boundary condition does nothing.

Boundary conditions are defined for each `meshField` as a `PtrList` of
`boundaryCondition`s. To correct the boundary conditions (i.e., to update the
ghost-cell values), `meshField`'s `correctBoundaryCondition()` function can be
called. This function updates ghost-cell values for each `meshLevel` and each
`meshDirection`. The implementation of the update is done in the `meshLevel`
class, such that one could also call the `correctBoundaryConditions()` function
on a specific level.

## Linear systems

When discretizing an equation in the FVM, one generally obtains a linear system
of equations of the form $Ax=b$, with $A$ the stencil matrix, $x$ the unknown
variable, and $b$ the right-hand side source. Such a linear system can be
defined in Briscola using the `linearSystem` class. This class is templated and
takes three template arguments: the stencil type (for $A$), the field type (for
$x$ and $b$) and the mesh type (for $A$, $x$ and $b$). The stencil type defines
a unique stencil size and shape which is used throughout the stencil matrix $A$.
Briscola's different stencil types are defined in

`src/briscolaCore/primitives/StencilSpace`.

The default `stencil` class (see `stencil.H`) has seven components: the center,
left, right, bottom, top, aft and fore components. For a given mesh cell, these
components correspond to the cell itself and its six direct neighbors across
cell faces. For each of these components, a scalar coefficient is defined. Also
the `diagStencil` class is defined, which trivially has just one diagonal
component. In principle, Briscola is designed to handle any sort of stencil,
possibly containing many more coefficients, but this is not yet implemented. The
`linearSystem` class can, given a solution field $x$, compute the residual $r
\equiv b-Ax$ as well as the evaluation $(Ax-b)/V$ with $V$ the cell volumes. The
division by $V$ is included because the discretizations contained in $A$ are
computed as volume integrals. In fact, whenever the right-hand side of the
linear system is modified, for example in:
```
linearSystem<stencil,vector,colocated> USys("sys(U)", U);
USys = im::ddt(U);
USys += source;
```
where `source` is an explicit source term, then 'under the hood', the addition
of `source` to `USys` is done as `b -= cv*source` with `cv` the cell volumes
(the subtraction is needed because `b` is in the right-hand side of the linear
system). See the `operator+=(...)` functions in
`src/briscolaFiniteVolume/linearSystems/linearSystem/linearSystem.C`

To obtain a solution to the linear system, it can be passed to one of the
several solvers available in Briscola. These are described in more detail in
[Chapter 7](./7_solvers.md).

## Schemes

Similarly to OpenFOAM, Briscola has several different generic classes to define
discretization schemes for spatial and temporal operators. These are defined in

`src/briscolaFiniteVolume/schemes`

and include discretizations for temporal derivatives (`ddtSchemes`), for
gradients (`gradientSchemes`), and many more. Discrete operator functions are
categorized as either:
* Implicit: see files included in `imSchemes.H`. These functions are part of the
`im` namespace (e.g., `im::div(...)`), which is analogous to the `fvm` namespace
in OpenFOAM. An implicit scheme returns a linear system with a stencil matrix
and a right-hand side source, such that a discrete linear system can be formed
by summing up several implicit (and optionally also explicit) discrete
operations (e.g., `im::ddt(u) + im::div(phi,u) = ...`). This is functionally
similar to how an `fvMatrix` can be constructed in OpenFOAM.
* Explicit: see files included in `exSchemes.H`. These functions are part of the
`ex` namespace (e.g., `ex::laplacian(...)`), which is analogous to the `fvc`
namespace in OpenFOAM. An explicit scheme returns a meshField which is equal to
the explicit evaluation of the chosen discrete operator.

In Briscola we've opted for the use of the `im` and `ex` namespaces as these
names are more meaningful and intuitive.

In addition to standard spatial and temporal operations, Briscola has schemes
for the following:
* Limiter schemes
* Prolongation and restriction schemes: these are needed for the geometric
multigrid solver, and define how fields can be prolonged or restricted to finer
and coarser grids respectively.
* Runge-Kutta schemes: For solvers which use Runge-Kutta-type schemes for
temporal discretization, the `RungeKuttaScheme` class and its derived
specializations define the number of intermediate stages and associated
coefficients in the form of a Butcher tableau (see Komen, E.M.J., et al.
"Analysis of the numerical dissipation rate of different Runge–Kutta and
velocity interpolation methods in an unstructured collocated finite volume
method in OpenFOAM." Computer Physics Communications 253 (2020): 107145.)

## Finite volume flow solvers

Briscola has several flow solvers which are based on the finite volume library.
These solver applications can be found in `applications/solvers`. They are
briefly described here:
* `briscolaColocated`: Colocated single-phase flow solver using Runge-Kutta time
discretization.
* `briscolaColocatedCNAB`: Colocated single-phase flow solver using
Crank-Nicolson discretization for the viscous term and Adamsh-Bashforth
discretization for the convective term.
* `briscolaColocatedTwoPhase`: Colocated two-phase flow solver using Runge-Kutta
time discretization and the Volume-of-Fluid (VOF) method for the two-phase
interface. More details on two-phase solvers are given in
[Chapter 5](./5_twoPhase.md).
* `briscolaColocatedTwoPhaseCNAB`: Colocated two-phase flow solver using
Crank-Nicolson discretization for the viscous term and Adamsh-Bashforth
discretization for the convective term.
* `briscolaLaplacian`: Colocated solver for a heat conduction equation of the
form $\partial T/\partial t=\nabla \cdot  (\alpha \nabla T) + S$, with temperature
$T$, diffusivity $\alpha$ and heat source $S$. The solver may also be used for
any other equation with the same form.
* `briscolaStaggered`: Staggered single-phase flow solver using Runge-Kutta time
discretization.
* `briscolaStaggeredCNAB`: Staggered single-phase flow solver using
Crank-Nicolson discretization for the viscous term and Adamsh-Bashforth
discretization for the convective term.
* `briscolaStaggeredTwoPhase`: Staggered two-phase flow solver using Runge-Kutta
time discretization.
* `briscolaStaggeredTwoPhaseCNAB`: Staggered two-phase flow solver using
Crank-Nicolson discretization for the viscous term and Adamsh-Bashforth
discretization for the convective term.
* `briscolaVofAdvection`: Colocated VOF advection solver without coupling to a
Navier-Stokes equation. Mainly used for verification test cases of VOF advection
under pre-defined velocity fields.

[Back to the table of contents](./0_start.md)
or [Next chapter: Two-phase solvers](./5_twoPhase.md)
