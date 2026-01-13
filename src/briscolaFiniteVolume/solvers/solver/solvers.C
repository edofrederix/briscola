#include "solvers.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolverBase(diagStencil,scalar,colocated);
makeSolverBase(diagStencil,scalar,staggered);
makeSolverBase(diagStencil,vector,colocated);
makeSolverBase(diagStencil,vector,staggered);

makeSolverBase(stencil,scalar,colocated);
makeSolverBase(stencil,scalar,staggered);
makeSolverBase(stencil,vector,colocated);
makeSolverBase(stencil,vector,staggered);

}

}

}
