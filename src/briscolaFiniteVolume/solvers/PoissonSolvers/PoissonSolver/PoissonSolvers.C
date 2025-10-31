#include "PoissonSolvers.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makePoissonSolverBase(stencil,scalar,colocated);
makePoissonSolverBase(stencil,scalar,staggered);
makePoissonSolverBase(stencil,vector,colocated);
makePoissonSolverBase(stencil,vector,staggered);

}

}

}
