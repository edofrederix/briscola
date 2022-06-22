#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
laplacianScheme<Type,MeshType>::laplacianScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type, class MeshType>
laplacianScheme<Type,MeshType>::laplacianScheme
(
    const laplacianScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
laplacianScheme<Type,MeshType>::~laplacianScheme()
{}

template<class Type, class MeshType>
autoPtr<laplacianScheme<Type,MeshType>>
laplacianScheme<Type,MeshType>::New
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

    return autoPtr<laplacianScheme<Type,MeshType>>(cstrIter()(dict, fvMsh));
}

}

}

}
