#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "APLU.H"
#include "Eigen.H"

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

makeDirectSolverType(Eigen,symmStencil,scalar,colocated);
makeDirectSolverType(Eigen,symmStencil,scalar,staggered);
makeDirectSolverType(Eigen,symmStencil,vector,colocated);
makeDirectSolverType(Eigen,symmStencil,vector,staggered);

makeDirectSolverType(Eigen,stencil,scalar,colocated);
makeDirectSolverType(Eigen,stencil,scalar,staggered);
makeDirectSolverType(Eigen,stencil,vector,colocated);
makeDirectSolverType(Eigen,stencil,vector,staggered);

}

}

}
