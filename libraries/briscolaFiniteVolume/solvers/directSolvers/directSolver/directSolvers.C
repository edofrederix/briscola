#include "directSolvers.H"
#include "addToRunTimeSelectionTable.H"

#include "APLU.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeDirectSolver(diagStencil,scalar,colocated);
makeDirectSolver(diagStencil,scalar,staggered);
makeDirectSolver(diagStencil,vector,colocated);
makeDirectSolver(diagStencil,vector,staggered);

makeDirectSolver(stencil,scalar,colocated);
makeDirectSolver(stencil,scalar,staggered);
makeDirectSolver(stencil,vector,colocated);
makeDirectSolver(stencil,vector,staggered);


makeDirectSolverType(APLU,diagStencil,scalar,colocated);
makeDirectSolverType(APLU,diagStencil,scalar,staggered);
makeDirectSolverType(APLU,diagStencil,vector,colocated);
makeDirectSolverType(APLU,diagStencil,vector,staggered);

makeDirectSolverType(APLU,stencil,scalar,colocated);
makeDirectSolverType(APLU,stencil,scalar,staggered);
makeDirectSolverType(APLU,stencil,vector,colocated);
makeDirectSolverType(APLU,stencil,vector,staggered);

}

}

}
