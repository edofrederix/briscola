#include "EigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(EigenSolver, 0);
defineRunTimeSelectionTable(EigenSolver, dictionary);

EigenSolver::EigenSolver()
{}

EigenSolver::EigenSolver(const EigenSolver& s)
{}

EigenSolver::~EigenSolver()
{}

autoPtr<EigenSolver> EigenSolver::New(const word solverType)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Eigen solver type " << solverType
            << ". Valid Eigen solvers are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<EigenSolver>(cstrIter()());
}

}

}

}
