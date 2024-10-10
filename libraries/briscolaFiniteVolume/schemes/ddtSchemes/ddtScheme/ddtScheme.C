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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "ddtSchemes", schemeName);

    word ddtSchemeType;
    is >> ddtSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(ddtSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown ddt scheme scheme "
            << ddtSchemeType << nl << nl
            << "Valid ddt schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ddtScheme<Type,MeshType>>(cstrIter()(fvMsh, is));
}

}

}

}
