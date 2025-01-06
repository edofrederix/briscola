#include "EigenSolver.H"

#include "SparseLUEigenSolver.H"
#include "PartialPivLUEigenSolver.H"
#include "BiCGSTABEigenSolver.H"

// The following solvers are only compiled when their macro variables are set.
// This is conditionally done in Make/options, based on the existence of certain
// environment variables.

#ifdef SUPERLU
#include "SuperLUEigenSolver.H"
#endif
