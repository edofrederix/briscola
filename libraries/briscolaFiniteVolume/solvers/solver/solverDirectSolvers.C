#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "APLU.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeDirectSolverType(APLU,diagStencil,scalar,colocated);
makeDirectSolverType(APLU,diagStencil,scalar,staggered);
makeDirectSolverType(APLU,diagStencil,vector,colocated);
makeDirectSolverType(APLU,diagStencil,vector,staggered);

makeDirectSolverType(APLU,symmStencil,scalar,colocated);
makeDirectSolverType(APLU,symmStencil,scalar,staggered);
makeDirectSolverType(APLU,symmStencil,vector,colocated);
makeDirectSolverType(APLU,symmStencil,vector,staggered);

makeDirectSolverType(APLU,stencil,scalar,colocated);
makeDirectSolverType(APLU,stencil,scalar,staggered);
makeDirectSolverType(APLU,stencil,vector,colocated);
makeDirectSolverType(APLU,stencil,vector,staggered);

}

}

}
