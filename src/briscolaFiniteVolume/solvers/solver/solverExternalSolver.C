#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::externalSolver::externalSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    dict_(dict),
    fvMsh_(fvMsh),
    l_(l),
    prepared_(false)
{}

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::externalSolver::~externalSolver()
{}

template<class SType, class Type, class MeshType>
autoPtr<typename solver<SType,Type,MeshType>::externalSolver>
solver<SType,Type,MeshType>::externalSolver::New
(
    const word& name,
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown external solver " << name << nl << nl
            << "Valid external solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solver<SType,Type,MeshType>::externalSolver>
    (
        cstrIter()(dict, fvMsh, l)
    );
}

template<class SType, class Type, class MeshType>
autoPtr<typename solver<SType,Type,MeshType>::externalSolver>
solver<SType,Type,MeshType>::externalSolver::New
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
{
    return solver<SType,Type,MeshType>::externalSolver::New
    (
        dict.lookup("type"),
        dict,
        fvMsh,
        l
    );
}

}

}

}
