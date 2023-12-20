#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
laplacianScheme<SType,Type,MeshType>::laplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class SType, class Type, class MeshType>
laplacianScheme<SType,Type,MeshType>::laplacianScheme
(
    const laplacianScheme<SType,Type,MeshType>& s
)
:
    scheme(s)
{}

template<class SType, class Type, class MeshType>
laplacianScheme<SType,Type,MeshType>::~laplacianScheme()
{}

template<class SType, class Type, class MeshType>
autoPtr<laplacianScheme<SType,Type,MeshType>>
laplacianScheme<SType,Type,MeshType>::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("laplacianSchemes").subDict(name)
    );

    const word laplacianSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(laplacianSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown laplacian scheme "
            << laplacianSchemeType << nl << nl
            << "Valid laplacian schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<laplacianScheme<SType,Type,MeshType>>
    (
        cstrIter()(dict, fvMsh)
    );
}

}

}

}
