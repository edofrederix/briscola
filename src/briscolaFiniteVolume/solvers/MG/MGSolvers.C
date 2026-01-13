#include "solvers.H"
#include "MG.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolverType(MG,diagStencil,scalar,colocated);
makeSolverType(MG,diagStencil,scalar,staggered);
makeSolverType(MG,diagStencil,vector,colocated);
makeSolverType(MG,diagStencil,vector,staggered);

makeSolverType(MG,stencil,scalar,colocated);
makeSolverType(MG,stencil,scalar,staggered);
makeSolverType(MG,stencil,vector,colocated);
makeSolverType(MG,stencil,vector,staggered);

}

}

}
