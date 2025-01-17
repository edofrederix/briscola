#include "EigenSolverBase.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(EigenSolverBase, 0);
defineRunTimeSelectionTable(EigenSolverBase, dictionary);

EigenSolverBase::EigenSolverBase(const dictionary& dict)
:
    dict_(dict)
{}

EigenSolverBase::EigenSolverBase(const EigenSolverBase& s)
:
    dict_(s.dict_)
{}

EigenSolverBase::~EigenSolverBase()
{}

autoPtr<EigenSolverBase> EigenSolverBase::New
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
            << "Unknown Eigen solver type " << solverType
            << ". Valid Eigen solvers are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<EigenSolverBase>(cstrIter()(dict));
}

void EigenSolverBase::prepare(const RowMajorMatrixType& A)
{
    NotImplemented;
}

void EigenSolverBase::prepare(const ColMajorMatrixType& A)
{
    NotImplemented;
}

}

}

}
