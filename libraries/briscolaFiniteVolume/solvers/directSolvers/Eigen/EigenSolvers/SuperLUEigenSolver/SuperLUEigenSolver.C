#include "SuperLUEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(SuperLUEigenSolver, 0);
addToRunTimeSelectionTable(EigenSolverBase, SuperLUEigenSolver, dictionary);

SuperLUEigenSolver::SuperLUEigenSolver(const dictionary& dict)
:
    EigenSolver<SuperLUSolver, ::Eigen::ColMajor>(dict)
{}

SuperLUEigenSolver::SuperLUEigenSolver(const SuperLUEigenSolver& s)
:
    EigenSolver<SuperLUSolver, ::Eigen::ColMajor>(s)
{}

SuperLUEigenSolver::~SuperLUEigenSolver()
{}

}

}

}
