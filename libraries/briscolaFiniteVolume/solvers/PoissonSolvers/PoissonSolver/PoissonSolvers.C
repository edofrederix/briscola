#include "PoissonSolvers.H"
#include "addToRunTimeSelectionTable.H"

#include "defaultPoissonSolver.H"
#include "FFTPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makePoissonSolver(stencil,scalar,colocated);
makePoissonSolver(stencil,scalar,staggered);
makePoissonSolver(stencil,vector,colocated);
makePoissonSolver(stencil,vector,staggered);


makePoissonSolverType(default,stencil,scalar,colocated);
makePoissonSolverType(default,stencil,scalar,staggered);
makePoissonSolverType(default,stencil,vector,colocated);
makePoissonSolverType(default,stencil,vector,staggered);

makePoissonSolverType(FFT,stencil,scalar,colocated);
makePoissonSolverType(FFT,stencil,scalar,staggered);
makePoissonSolverType(FFT,stencil,vector,colocated);
makePoissonSolverType(FFT,stencil,vector,staggered);

}

}

}
