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

void PardisoEigenSolver::prepare
(
    const EigenMatrixType& A,
    const scalar
)
{
    solverPtr_.reset(new EigenSolverType(A));
    solverPtr_->compute(A);
}

void PardisoEigenSolver::solve(EigenRhsType& x, const EigenRhsType& b) const
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
