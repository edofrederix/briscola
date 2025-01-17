#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::smoother::smoother
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    dict_(dict),
    fvMsh_(fvMsh)
{}

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::smoother::~smoother()
{}

template<class SType, class Type, class MeshType>
autoPtr<typename solver<SType,Type,MeshType>::smoother>
solver<SType,Type,MeshType>::smoother::New
(
    const word& name,
    const dictionary& dict,
    const fvMesh& fvMsh
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown smoother " << name << nl << nl
            << "Valid smoothers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solver<SType,Type,MeshType>::smoother>
    (
        cstrIter()(dict, fvMsh)
    );
}

}

}

}
