#include "PCLUPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(PCLUPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, PCLUPETScSolver, dictionary);

PCLUPETScSolver::PCLUPETScSolver(const dictionary& dict)
:
    PETScExternalSolver(dict)
{}

PCLUPETScSolver::PCLUPETScSolver(const PCLUPETScSolver& s)
:
    PETScExternalSolver(s)
{}

PCLUPETScSolver::~PCLUPETScSolver()
{}

void PCLUPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScExternalSolver::prepare(mat, comm);

    Mat& f = factor_();

    MatGetFactor(mat, MATSOLVERPETSC, MAT_FACTOR_LU, &f);

    int m, n;
    MatGetLocalSize(mat, &m, &n);

    IS row_perm, col_perm;

    // Create identity permutations (no reordering)

    ISCreateStride(comm, m, 0, 1, &row_perm);
    ISCreateStride(comm, m, 0, 1, &col_perm);

    MatLUFactorSymbolic(f, mat, row_perm, col_perm, NULL);
    MatLUFactorNumeric(f, mat, NULL);
}

}

}

}
