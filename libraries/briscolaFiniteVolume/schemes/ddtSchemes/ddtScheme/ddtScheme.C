#include "ddtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
ddtScheme<Type,MeshType>::ddtScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type, class MeshType>
ddtScheme<Type,MeshType>::ddtScheme
(
    const ddtScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
ddtScheme<Type,MeshType>::~ddtScheme()
{}

template<class Type, class MeshType>
autoPtr<ddtScheme<Type,MeshType>> ddtScheme<Type,MeshType>::New
(
    const word schemeName,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("ddtSchemes").subDict(schemeName)
    );

    const word ddtSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(ddtSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown ddt scheme scheme "
            << ddtSchemeType << nl << nl
            << "Valid ddt schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ddtScheme<Type,MeshType>>(cstrIter()(dict, fvMsh));
}

}

}

}
