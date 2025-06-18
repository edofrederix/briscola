# Chapter 8: Current limitations

[Back to the table of contents](./0_start.md)

Briscola is in development, and currently has some limitations which are listed
here. Briscola's development team is small and its resources are limited. Thus,
contributions from the community towards resolving these limitations are
encouraged and can be proposed via Github.

* *Two-phase solver:* While most normal and interface reconstruction algorithms
  are capable of handling general hexahedron cells, the split advection solver
  is currently not able to handle unstructured meshes.
* *FFT solver:* The split FFT solver only works well for mild density ratios.
  Improving this requires fundamental improvements to the FFT algorithm, or to
  the pressure extrapolation method.
* *Staggered solvers:* Staggered solvers only work on rectilinear meshes that
  are aligned with the coordinate system. Following the discretizations proposed
  in Wesseling, P. Principles of computational fluid dynamics. Vol. 29. Springer
  Science & Business Media, 2009., this may be generalized to arbitrary grids
  but that is yet to be implemented.
* *GPU computing:* There is currently no GPU computing support
* *Ghost cells:* The number of ghost cell layers can be generalized, in order to
  flexibly support wider stencils.
* *Two-phase solver:* The surface tension approach is explicit, and can impose a
  restriction on the time step size. This is currently not checked, and may lead
  to instability.
* ...

Feel free to report any additional issues via Github or by contacting the
authors directly.

[Back to the table of contents](./0_start.md)
