#include "solvers.H"
#include "Eigen.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeExternalSolverType(Eigen,stencil,scalar,colocated);
makeExternalSolverType(Eigen,stencil,scalar,staggered);
makeExternalSolverType(Eigen,stencil,vector,colocated);
makeExternalSolverType(Eigen,stencil,vector,staggered);

}

}

}
