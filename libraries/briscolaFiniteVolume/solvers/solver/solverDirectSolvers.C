#include "solvers.H"
#include "addToRunTimeSelectionTable.H"

#include "Eigen.H"
#include "PETSc.H"

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

makeDirectSolverType(PETSc,symmStencil,scalar,colocated);
makeDirectSolverType(PETSc,symmStencil,scalar,staggered);
makeDirectSolverType(PETSc,symmStencil,vector,colocated);
makeDirectSolverType(PETSc,symmStencil,vector,staggered);

makeDirectSolverType(PETSc,stencil,scalar,colocated);
makeDirectSolverType(PETSc,stencil,scalar,staggered);
makeDirectSolverType(PETSc,stencil,vector,colocated);
makeDirectSolverType(PETSc,stencil,vector,staggered);

}

}

}
