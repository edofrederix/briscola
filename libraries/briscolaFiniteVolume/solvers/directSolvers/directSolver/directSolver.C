#include "directSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
directSolver<SType,Type,MeshType>::directSolver
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    dict_(dict),
    fvMsh_(fvMsh),
    l_(l)
{}

template<class SType, class Type, class MeshType>
directSolver<SType,Type,MeshType>::~directSolver()
{}

template<class SType, class Type, class MeshType>
autoPtr<directSolver<SType,Type,MeshType>>
directSolver<SType,Type,MeshType>::New
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
{
    const word solverType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown direct solver "
            << solverType << nl << nl
            << "Valid direct solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<directSolver<SType,Type,MeshType>>
    (
        cstrIter()(dict, fvMsh, l)
    );
}

}

}

}
