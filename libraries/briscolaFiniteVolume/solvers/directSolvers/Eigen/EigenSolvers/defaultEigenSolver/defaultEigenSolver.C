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

void defaultEigenSolver::prepare(const matrixType& A)
{
    if (isPositiveDefinite(A))
    {
        gSolverPtr_.clear();
        sSolverPtr_.reset(new sSolverType(A));
        sSolverPtr_->compute(A);
    }
    else
    {
        sSolverPtr_.clear();
        gSolverPtr_.reset(new gSolverType(A));
        gSolverPtr_->compute(A);
    }
}

void defaultEigenSolver::solve(rhsType& x, const rhsType& b) const
{
    if (gSolverPtr_.valid())
    {
        x = gSolverPtr_->solve(b);
    }
    else if (sSolverPtr_.valid())
    {
        x = sSolverPtr_->solve(b);
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
