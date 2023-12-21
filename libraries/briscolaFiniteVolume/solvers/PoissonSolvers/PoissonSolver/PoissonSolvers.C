#include "PoissonSolvers.H"
#include "addToRunTimeSelectionTable.H"

#include "defaultPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makePoissonSolver(stencil,scalar,colocated);
makePoissonSolver(stencil,scalar,staggered);
makePoissonSolver(stencil,vector,colocated);
makePoissonSolver(stencil,vector,staggered);

makePoissonSolver(symmStencil,scalar,colocated);
makePoissonSolver(symmStencil,scalar,staggered);
makePoissonSolver(symmStencil,vector,colocated);
makePoissonSolver(symmStencil,vector,staggered);

makePoissonSolverType(default,stencil,scalar,colocated);
makePoissonSolverType(default,stencil,vector,colocated);
makePoissonSolverType(default,symmStencil,scalar,colocated);
makePoissonSolverType(default,symmStencil,vector,colocated);

makePoissonSolverType(default,stencil,scalar,staggered);
makePoissonSolverType(default,stencil,vector,staggered);
makePoissonSolverType(default,symmStencil,scalar,staggered);
makePoissonSolverType(default,symmStencil,vector,staggered);

}

}

}
