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

makeExternalSolverType(Eigen,stencil,scalar,colocated);
makeExternalSolverType(Eigen,stencil,scalar,staggered);
makeExternalSolverType(Eigen,stencil,vector,colocated);
makeExternalSolverType(Eigen,stencil,vector,staggered);

makeExternalSolverType(PETSc,stencil,scalar,colocated);
makeExternalSolverType(PETSc,stencil,scalar,staggered);
makeExternalSolverType(PETSc,stencil,vector,colocated);
makeExternalSolverType(PETSc,stencil,vector,staggered);

}

}

}
