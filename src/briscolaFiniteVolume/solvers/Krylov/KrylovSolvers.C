#ifdef PETSC

#include "solvers.H"
#include "Krylov.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolverType(Krylov,diagStencil,scalar,colocated);
makeSolverType(Krylov,diagStencil,scalar,staggered);
makeSolverType(Krylov,diagStencil,vector,colocated);
makeSolverType(Krylov,diagStencil,vector,staggered);

makeSolverType(Krylov,stencil,scalar,colocated);
makeSolverType(Krylov,stencil,scalar,staggered);
makeSolverType(Krylov,stencil,vector,colocated);
makeSolverType(Krylov,stencil,vector,staggered);

}

}

}

#endif
