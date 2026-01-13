#include "PoissonSolvers.H"
#include "defaultPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makePoissonSolverType(default,stencil,scalar,colocated);
makePoissonSolverType(default,stencil,vector,colocated);

makePoissonSolverType(default,stencil,scalar,staggered);
makePoissonSolverType(default,stencil,vector,staggered);

}

}

}
