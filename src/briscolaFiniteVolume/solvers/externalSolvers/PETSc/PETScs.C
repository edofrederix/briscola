#ifdef PETSC

#include "solvers.H"
#include "PETSc.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeExternalSolverType(PETSc,stencil,scalar,colocated);
makeExternalSolverType(PETSc,stencil,scalar,staggered);
makeExternalSolverType(PETSc,stencil,vector,colocated);
makeExternalSolverType(PETSc,stencil,vector,staggered);

}

}

}

#endif
