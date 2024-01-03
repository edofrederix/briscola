#include "PardisoEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(PardisoEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolver, PardisoEigenSolver, dictionary);

PardisoEigenSolver::PardisoEigenSolver()
:
    EigenSolver()
{
    ::Eigen::initParallel();
}

PardisoEigenSolver::PardisoEigenSolver(const PardisoEigenSolver& s)
:
    EigenSolver(s)
{
    ::Eigen::initParallel();
}

PardisoEigenSolver::~PardisoEigenSolver()
{}

void PardisoEigenSolver::prepare(const matrixType& A)
{
    const bool symm = isSymmetric(A);

    if (symm)
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

void PardisoEigenSolver::solve(rhsType& x, const rhsType& b) const
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
