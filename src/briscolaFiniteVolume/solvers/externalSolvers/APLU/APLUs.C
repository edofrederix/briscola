#include "solvers.H"
#include "APLU.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeExternalSolverType(APLU,stencil,scalar,colocated);
makeExternalSolverType(APLU,stencil,scalar,staggered);
makeExternalSolverType(APLU,stencil,vector,colocated);
makeExternalSolverType(APLU,stencil,vector,staggered);

}

}

}
