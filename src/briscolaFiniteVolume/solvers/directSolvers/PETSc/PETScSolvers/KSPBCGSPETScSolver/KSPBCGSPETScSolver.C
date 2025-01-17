#include "KSPBCGSPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(KSPBCGSPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, KSPBCGSPETScSolver, dictionary);

KSPBCGSPETScSolver::KSPBCGSPETScSolver(const dictionary& dict)
:
    PETScKSPSolver(dict)
{}

KSPBCGSPETScSolver::KSPBCGSPETScSolver(const KSPBCGSPETScSolver& s)
:
    PETScKSPSolver(s)
{}

KSPBCGSPETScSolver::~KSPBCGSPETScSolver()
{}

void KSPBCGSPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScKSPSolver::prepare(mat, comm);
    KSPSetType(solverPtr_(), KSPBCGS);
}

}

}

}
