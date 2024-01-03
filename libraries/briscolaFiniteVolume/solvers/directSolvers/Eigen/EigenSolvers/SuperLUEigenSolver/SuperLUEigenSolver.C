#include "SuperLUEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SuperLUEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, SuperLUEigenSolver, dictionary);

SuperLUEigenSolver::SuperLUEigenSolver()
:
    EigenSolver()
{}

SuperLUEigenSolver::SuperLUEigenSolver(const SuperLUEigenSolver& s)
:
    EigenSolver(s)
{}

SuperLUEigenSolver::~SuperLUEigenSolver()
{}

void SuperLUEigenSolver::prepare(const matrixType& A)
{
    solverPtr_.reset(new solverType(A));
    solverPtr_->compute(A);
}

void SuperLUEigenSolver::solve(rhsType& x, const rhsType& b) const
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
