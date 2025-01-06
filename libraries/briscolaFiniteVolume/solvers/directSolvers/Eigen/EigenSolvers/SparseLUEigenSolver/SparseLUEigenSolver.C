#include "SparseLUEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SparseLUEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolverBase, SparseLUEigenSolver, dictionary);

SparseLUEigenSolver::SparseLUEigenSolver(const dictionary& dict)
:
    EigenSolver<SparseLUSolver, ::Eigen::ColMajor>(dict)
{}

SparseLUEigenSolver::SparseLUEigenSolver(const SparseLUEigenSolver& s)
:
    EigenSolver<SparseLUSolver, ::Eigen::ColMajor>(s)
{}

SparseLUEigenSolver::~SparseLUEigenSolver()
{}

}

}

}
