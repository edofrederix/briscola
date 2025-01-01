#include "BiCGSTABEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(BiCGSTABEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, BiCGSTABEigenSolver, dictionary);

BiCGSTABEigenSolver::BiCGSTABEigenSolver()
:
    EigenSolver()
{}

BiCGSTABEigenSolver::BiCGSTABEigenSolver(const BiCGSTABEigenSolver& s)
:
    EigenSolver(s)
{}

BiCGSTABEigenSolver::~BiCGSTABEigenSolver()
{}

void BiCGSTABEigenSolver::prepare
(
    const EigenMatrixType& A,
    const scalar tolerance
)
{
    solverPtr_.reset(new EigenSolverType(A));
    solverPtr_->setTolerance(tolerance);
    solverPtr_->compute(A);
}

void BiCGSTABEigenSolver::solve(EigenRhsType& x, const EigenRhsType& b) const
{
    if (solverPtr_.valid())
    {
        x = solverPtr_->solveWithGuess(b, x);
    }
    else
    {
        FatalErrorInFunction
            << "Solver not prepared" << endl << abort(FatalError);
    }
}

}

}

}
