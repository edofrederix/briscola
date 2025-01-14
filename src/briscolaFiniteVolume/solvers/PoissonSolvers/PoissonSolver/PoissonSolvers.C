#include "PoissonSolvers.H"
#include "addToRunTimeSelectionTable.H"

#include "splitPoissonSolver.H"
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

makePoissonSolverType(default,stencil,scalar,colocated);
makePoissonSolverType(default,stencil,vector,colocated);

makePoissonSolverType(default,stencil,scalar,staggered);
makePoissonSolverType(default,stencil,vector,staggered);

makePoissonSolverType(split,stencil,scalar,colocated);
makePoissonSolverType(split,stencil,vector,colocated);

makePoissonSolverType(split,stencil,scalar,staggered);
makePoissonSolverType(split,stencil,vector,staggered);

}

}

}
