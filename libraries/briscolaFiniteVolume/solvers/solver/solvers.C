#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "MGSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolver(diagStencil,scalar,colocated);
makeSolver(diagStencil,scalar,staggered);
makeSolver(diagStencil,vector,colocated);
makeSolver(diagStencil,vector,staggered);

makeSolver(stencil,scalar,colocated);
makeSolver(stencil,scalar,staggered);
makeSolver(stencil,vector,colocated);
makeSolver(stencil,vector,staggered);


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
