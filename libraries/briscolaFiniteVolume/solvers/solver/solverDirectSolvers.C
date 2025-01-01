#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "Eigen.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

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
