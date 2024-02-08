#include "defaultEigenSolver.H"

// The following solvers are only compiled when their macro variables are set.
// This is conditionally done in Make/options, based on the existance of certain
// environment variables.

#ifdef MKL
#include "PardisoEigenSolver.H"
#endif

#ifdef SUITESPARSE
#include "UmfPackEigenSolver.H"
#endif

#ifdef SUPERLU
#include "SuperLUEigenSolver.H"
#endif
