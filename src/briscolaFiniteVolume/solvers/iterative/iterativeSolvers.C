#include "solvers.H"
#include "iterative.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolverType(iterative,diagStencil,scalar,colocated);
makeSolverType(iterative,diagStencil,scalar,staggered);
makeSolverType(iterative,diagStencil,vector,colocated);
makeSolverType(iterative,diagStencil,vector,staggered);

makeSolverType(iterative,stencil,scalar,colocated);
makeSolverType(iterative,stencil,scalar,staggered);
makeSolverType(iterative,stencil,vector,colocated);
makeSolverType(iterative,stencil,vector,staggered);

}

}

}
