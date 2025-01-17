#include "SuperLUDistPETScSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SuperLUDistPETScSolver, 0);
addToRunTimeSelectionTable(PETScSolverBase, SuperLUDistPETScSolver, dictionary);

SuperLUDistPETScSolver::SuperLUDistPETScSolver(const dictionary& dict)
:
    PETScDirectSolver(dict)
{}

SuperLUDistPETScSolver::SuperLUDistPETScSolver(const SuperLUDistPETScSolver& s)
:
    PETScDirectSolver(s)
{}

SuperLUDistPETScSolver::~SuperLUDistPETScSolver()
{}

void SuperLUDistPETScSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    PETScDirectSolver::prepare(mat, comm);

    Mat& f = factor_();

    MatGetFactor(mat, MATSOLVERSUPERLU_DIST, MAT_FACTOR_LU, &f);

    MatLUFactorSymbolic(f, mat, NULL, NULL, NULL);
    MatLUFactorNumeric(f, mat, NULL);
}

}

}

}
