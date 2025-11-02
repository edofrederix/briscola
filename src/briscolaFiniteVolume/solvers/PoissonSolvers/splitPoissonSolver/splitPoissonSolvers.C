#include "PoissonSolvers.H"
#include "splitPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makePoissonSolverType(split,stencil,scalar,colocated);
makePoissonSolverType(split,stencil,vector,colocated);

makePoissonSolverType(split,stencil,scalar,staggered);
makePoissonSolverType(split,stencil,vector,staggered);

}

}

}
