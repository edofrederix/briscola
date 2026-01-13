#include "solvers.H"
#include "diagonal.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolverType(diagonal,diagStencil,scalar,colocated);
makeSolverType(diagonal,diagStencil,scalar,staggered);
makeSolverType(diagonal,diagStencil,vector,colocated);
makeSolverType(diagonal,diagStencil,vector,staggered);

makeSolverType(diagonal,stencil,scalar,colocated);
makeSolverType(diagonal,stencil,scalar,staggered);
makeSolverType(diagonal,stencil,vector,colocated);
makeSolverType(diagonal,stencil,vector,staggered);

}

}

}
