#include "KSPIBCGSPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(KSPIBCGSPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, KSPIBCGSPETScSolver, dictionary);

KSPIBCGSPETScSolver::KSPIBCGSPETScSolver(const dictionary& dict)
:
    PETScKSPSolver(dict)
{}

KSPIBCGSPETScSolver::KSPIBCGSPETScSolver(const KSPIBCGSPETScSolver& s)
:
    PETScKSPSolver(s)
{}

KSPIBCGSPETScSolver::~KSPIBCGSPETScSolver()
{}

void KSPIBCGSPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScKSPSolver::prepare(mat, comm);
    KSPSetType(solverPtr_(), KSPIBCGS);
}

}

}

}
