#include "partialPivLUEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(partialPivLUEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, partialPivLUEigenSolver, dictionary);

partialPivLUEigenSolver::partialPivLUEigenSolver()
:
    EigenSolver()
{}

partialPivLUEigenSolver::partialPivLUEigenSolver(const partialPivLUEigenSolver& s)
:
    EigenSolver(s)
{}

partialPivLUEigenSolver::~partialPivLUEigenSolver()
{}

void partialPivLUEigenSolver::prepare
(
    const EigenMatrixType& A,
    const scalar
)
{
    solverPtr_.reset(new EigenSolverType(A));
    solverPtr_->compute(A);
}

void partialPivLUEigenSolver::solve(EigenRhsType& x, const EigenRhsType& b)
const
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
