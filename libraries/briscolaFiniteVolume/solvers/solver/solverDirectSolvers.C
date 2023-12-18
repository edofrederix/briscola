#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "APLU.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

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
