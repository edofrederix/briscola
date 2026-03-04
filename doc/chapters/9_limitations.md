# Chapter 9: Current limitations

[Back to the table of contents](./0_start.md)
or [Previous chapter: Testing and validation](./8_testing.md)

Briscola is in development, and currently has some limitations which are listed
here. Briscola's development team is small and its resources are limited. Thus,
contributions from the community towards resolving these limitations are
encouraged and can be proposed via GitHub.

* *Two-phase solver:* While most normal and interface reconstruction algorithms
  can handle general hexahedron cells, the split advection solver is currently
  not able to handle unstructured meshes.
* *FFT solver:* For two-phase simulations, the split FFT solver only works well
  for mild density ratios. Improving this requires fundamental improvements to
  the solver's algorithm, or to the pressure extrapolation method.
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
* *Immersed boundaries*: Currently, only infinite cylinders and spheres (both
  normal and inverted) can be used to define immersed boundaries. Additional
  analytical shape definitions are yet to be implemented. Defining immersed
  boundaries based on CAD files is also currently not supported.
* When the `FOAM_SIGFPE` environment variable is set, sigFpe trapping is enabled
  and this may lead to problems in Briscola. The reason is that ghost cell
  values are sometimes not initialized, but they are treated in some operations
  anyway. As a result invalid values may be accessed. Usually this is harmless,
  but with `FOAM_SIGFPE` set the code will crash anyway. Consider disabling
  `FOAM_SIGFPE` (i.e., `unset FOAM_SIGFPE`), or handle special situations using
  the `sigFpeEnabled()` function, as is for example done in a some schemes and
  test applications.

Feel free to report any additional issues via GitHub or by contacting the
authors directly.

[Back to the table of contents](./0_start.md)
