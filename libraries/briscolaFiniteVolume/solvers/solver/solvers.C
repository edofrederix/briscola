#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "MG.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeSolver(symmStencil,scalar,colocated);
makeSolver(symmStencil,scalar,staggered);
makeSolver(symmStencil,vector,colocated);
makeSolver(symmStencil,vector,staggered);

makeSolver(stencil,scalar,colocated);
makeSolver(stencil,scalar,staggered);
makeSolver(stencil,vector,colocated);
makeSolver(stencil,vector,staggered);


makeSolverType(MG,symmStencil,scalar,colocated);
makeSolverType(MG,symmStencil,scalar,staggered);
makeSolverType(MG,symmStencil,vector,colocated);
makeSolverType(MG,symmStencil,vector,staggered);

makeSolverType(MG,stencil,scalar,colocated);
makeSolverType(MG,stencil,scalar,staggered);
makeSolverType(MG,stencil,vector,colocated);
makeSolverType(MG,stencil,vector,staggered);

}

}

}
