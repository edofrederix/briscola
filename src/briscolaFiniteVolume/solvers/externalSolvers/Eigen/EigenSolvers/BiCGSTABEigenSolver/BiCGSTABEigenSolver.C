#include "BiCGSTABEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(BiCGSTABEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolverBase, BiCGSTABEigenSolver, dictionary);

BiCGSTABEigenSolver::BiCGSTABEigenSolver(const dictionary& dict)
:
    EigenSolver<BiCGSTABSolver,::Eigen::RowMajor>(dict),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 1000)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", -1))
{}

BiCGSTABEigenSolver::BiCGSTABEigenSolver(const BiCGSTABEigenSolver& s)
:
    EigenSolver<BiCGSTABSolver,::Eigen::RowMajor>(s),
    maxIter_(s.maxIter_)
{}

BiCGSTABEigenSolver::~BiCGSTABEigenSolver()
{}

void BiCGSTABEigenSolver::solve(RhsType& x, const RhsType& b)
{
    if (solverPtr_.valid())
    {
        if (tolerance_ > 0)
            solverPtr_->setTolerance(tolerance_);

        solverPtr_->setMaxIterations(maxIter_);

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
