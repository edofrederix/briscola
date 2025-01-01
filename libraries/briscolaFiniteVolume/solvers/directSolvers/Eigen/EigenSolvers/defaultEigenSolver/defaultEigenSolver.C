#include "defaultEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(defaultEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, defaultEigenSolver, dictionary);

defaultEigenSolver::defaultEigenSolver()
:
    EigenSolver()
{}

defaultEigenSolver::defaultEigenSolver(const defaultEigenSolver& s)
:
    EigenSolver(s)
{}

defaultEigenSolver::~defaultEigenSolver()
{}

void defaultEigenSolver::prepare
(
    const EigenMatrixType& A,
    const scalar
)
{
    solverPtr_.reset(new EigenSolverType(A));
    solverPtr_->compute(A);
}

void defaultEigenSolver::solve(EigenRhsType& x, const EigenRhsType& b) const
{
    if (solverPtr_.valid())
    {
        x = solverPtr_->solve(b);
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
