#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "MG.H"
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

makeSolver(symmStencil,scalar,colocated);
makeSolver(symmStencil,scalar,staggered);
makeSolver(symmStencil,vector,colocated);
makeSolver(symmStencil,vector,staggered);

makeSolver(stencil,scalar,colocated);
makeSolver(stencil,scalar,staggered);
makeSolver(stencil,vector,colocated);
makeSolver(stencil,vector,staggered);


makeSolverType(MG,diagStencil,scalar,colocated);
makeSolverType(MG,diagStencil,scalar,staggered);
makeSolverType(MG,diagStencil,vector,colocated);
makeSolverType(MG,diagStencil,vector,staggered);

makeSolverType(MG,symmStencil,scalar,colocated);
makeSolverType(MG,symmStencil,scalar,staggered);
makeSolverType(MG,symmStencil,vector,colocated);
makeSolverType(MG,symmStencil,vector,staggered);

makeSolverType(MG,stencil,scalar,colocated);
makeSolverType(MG,stencil,scalar,staggered);
makeSolverType(MG,stencil,vector,colocated);
makeSolverType(MG,stencil,vector,staggered);


makeSolverType(diagonal,diagStencil,scalar,colocated);
makeSolverType(diagonal,diagStencil,scalar,staggered);
makeSolverType(diagonal,diagStencil,vector,colocated);
makeSolverType(diagonal,diagStencil,vector,staggered);

makeSolverType(diagonal,symmStencil,scalar,colocated);
makeSolverType(diagonal,symmStencil,scalar,staggered);
makeSolverType(diagonal,symmStencil,vector,colocated);
makeSolverType(diagonal,symmStencil,vector,staggered);

makeSolverType(diagonal,stencil,scalar,colocated);
makeSolverType(diagonal,stencil,scalar,staggered);
makeSolverType(diagonal,stencil,vector,colocated);
makeSolverType(diagonal,stencil,vector,staggered);

}

}

}
