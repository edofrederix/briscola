#include "KSPFGMRESPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(KSPFGMRESPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, KSPFGMRESPETScSolver, dictionary);

KSPFGMRESPETScSolver::KSPFGMRESPETScSolver(const dictionary& dict)
:
    PETScKSPSolver(dict)
{}

KSPFGMRESPETScSolver::KSPFGMRESPETScSolver(const KSPFGMRESPETScSolver& s)
:
    PETScKSPSolver(s)
{}

KSPFGMRESPETScSolver::~KSPFGMRESPETScSolver()
{}

void KSPFGMRESPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScKSPSolver::prepare(mat, comm);
    KSPSetType(solverPtr_(), KSPFGMRES);
}

}

}

}
