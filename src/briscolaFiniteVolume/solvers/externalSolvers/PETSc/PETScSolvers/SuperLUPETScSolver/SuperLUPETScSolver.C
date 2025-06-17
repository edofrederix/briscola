#include "SuperLUPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SuperLUPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, SuperLUPETScSolver, dictionary);

SuperLUPETScSolver::SuperLUPETScSolver(const dictionary& dict)
:
    PETScExternalSolver(dict)
{}

SuperLUPETScSolver::SuperLUPETScSolver(const SuperLUPETScSolver& s)
:
    PETScExternalSolver(s)
{}

SuperLUPETScSolver::~SuperLUPETScSolver()
{}

void SuperLUPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScExternalSolver::prepare(mat, comm);

    Mat& f = factor_();

    MatGetFactor(mat, MATSOLVERSUPERLU, MAT_FACTOR_LU, &f);

    MatLUFactorSymbolic(f, mat, NULL, NULL, NULL);
    MatLUFactorNumeric(f, mat, NULL);
}

}

}

}
