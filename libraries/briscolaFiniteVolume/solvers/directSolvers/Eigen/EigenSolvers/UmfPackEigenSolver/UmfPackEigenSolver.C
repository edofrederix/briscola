#include "UmfPackEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(UmfPackEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, UmfPackEigenSolver, dictionary);

UmfPackEigenSolver::UmfPackEigenSolver()
:
    EigenSolver()
{}

UmfPackEigenSolver::UmfPackEigenSolver(const UmfPackEigenSolver& s)
:
    EigenSolver(s)
{}

UmfPackEigenSolver::~UmfPackEigenSolver()
{}

void UmfPackEigenSolver::prepare
(
    const EigenMatrixType& A,
    const scalar
)
{
    solverPtr_.reset(new EigenSolverType(A));
    solverPtr_->compute(A);
}

void UmfPackEigenSolver::solve(EigenRhsType& x, const EigenRhsType& b) const
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
