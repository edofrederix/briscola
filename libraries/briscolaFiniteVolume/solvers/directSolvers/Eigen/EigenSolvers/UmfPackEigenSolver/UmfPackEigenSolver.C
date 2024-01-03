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

void UmfPackEigenSolver::prepare(const matrixType& A)
{
    solverPtr_.reset(new solverType(A));
    solverPtr_->compute(A);
}

void UmfPackEigenSolver::solve(rhsType& x, const rhsType& b) const
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
