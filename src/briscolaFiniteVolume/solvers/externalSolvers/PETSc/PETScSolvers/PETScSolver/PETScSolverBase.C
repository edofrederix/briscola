#include "PETScSolverBase.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(PETScSolverBase, 0);
defineRunTimeSelectionTable(PETScSolverBase, dictionary);

PETScSolverBase::PETScSolverBase(const dictionary& dict)
:
    dict_(dict)
{}

PETScSolverBase::PETScSolverBase(const PETScSolverBase& s)
:
    dict_(s.dict_)
{}

PETScSolverBase::~PETScSolverBase()
{}

autoPtr<PETScSolverBase> PETScSolverBase::New
(
    const word solverType,
    const dictionary& dict
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown PETSc solver type " << solverType
            << ". Valid PETSc solvers are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<PETScSolverBase>(cstrIter()(dict));
}

}

}

}
