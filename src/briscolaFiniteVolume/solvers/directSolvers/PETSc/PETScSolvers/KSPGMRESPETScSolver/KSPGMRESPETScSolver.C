#include "KSPGMRESPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(KSPGMRESPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, KSPGMRESPETScSolver, dictionary);

KSPGMRESPETScSolver::KSPGMRESPETScSolver(const dictionary& dict)
:
    PETScKSPSolver(dict)
{}

KSPGMRESPETScSolver::KSPGMRESPETScSolver(const KSPGMRESPETScSolver& s)
:
    PETScKSPSolver(s)
{}

KSPGMRESPETScSolver::~KSPGMRESPETScSolver()
{}

void KSPGMRESPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScKSPSolver::prepare(mat, comm);
    KSPSetType(solverPtr_(), KSPGMRES);
}

}

}

}
