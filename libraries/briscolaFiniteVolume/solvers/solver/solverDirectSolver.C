#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::directSolver::directSolver
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
solver<SType,Type,MeshType>::directSolver::~directSolver()
{}

template<class SType, class Type, class MeshType>
autoPtr<typename solver<SType,Type,MeshType>::directSolver>
solver<SType,Type,MeshType>::directSolver::New
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
            << "Unknown direct solver " << name << nl << nl
            << "Valid direct solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solver<SType,Type,MeshType>::directSolver>
    (
        cstrIter()(dict, fvMsh, l)
    );
}

template<class SType, class Type, class MeshType>
autoPtr<typename solver<SType,Type,MeshType>::directSolver>
solver<SType,Type,MeshType>::directSolver::New
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
{
    return solver<SType,Type,MeshType>::directSolver::New
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
