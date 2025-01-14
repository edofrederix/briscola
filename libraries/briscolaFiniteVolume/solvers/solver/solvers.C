#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "MG.H"
#include "Krylov.H"
#include "diagonal.H"

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


makeSolverType(Krylov,diagStencil,scalar,colocated);
makeSolverType(Krylov,diagStencil,scalar,staggered);
makeSolverType(Krylov,diagStencil,vector,colocated);
makeSolverType(Krylov,diagStencil,vector,staggered);

makeSolverType(Krylov,stencil,scalar,colocated);
makeSolverType(Krylov,stencil,scalar,staggered);
makeSolverType(Krylov,stencil,vector,colocated);
makeSolverType(Krylov,stencil,vector,staggered);


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
