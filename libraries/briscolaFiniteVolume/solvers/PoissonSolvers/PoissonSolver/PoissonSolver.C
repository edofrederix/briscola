#include "PoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
PoissonSolver<SType,Type,MeshType>::PoissonSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    dict_(dict),
    fvMsh_(fvMsh)
{}

template<class SType, class Type, class MeshType>
PoissonSolver<SType,Type,MeshType>::~PoissonSolver()
{}

template<class SType, class Type, class MeshType>
autoPtr<PoissonSolver<SType,Type,MeshType>> PoissonSolver<SType,Type,MeshType>::New
(
    const word PoissonSolverName,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.solverDict().subDict(PoissonSolverName)
    );

    const word PoissonSolverType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(PoissonSolverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find("default");
    }

    return autoPtr<PoissonSolver<SType,Type,MeshType>>
    (
        cstrIter()(PoissonSolverName, dict, fvMsh)
    );
}

}

}

}
