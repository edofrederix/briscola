#include "PartialPivLUEigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(PartialPivLUEigenSolver, 0);
addToRunTimeSelectionTable
(
    EigenSolverBase,
    PartialPivLUEigenSolver,
    dictionary
);

PartialPivLUEigenSolver::PartialPivLUEigenSolver(const dictionary& dict)
:
    EigenSolver<PartialPivLUSolver, ::Eigen::ColMajor>(dict)
{}

PartialPivLUEigenSolver::PartialPivLUEigenSolver
(
    const PartialPivLUEigenSolver& s
)
:
    EigenSolver<PartialPivLUSolver, ::Eigen::ColMajor>(s)
{}

PartialPivLUEigenSolver::~PartialPivLUEigenSolver()
{}

}

}

}
