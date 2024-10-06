#include "faceGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
faceGradientScheme<Type,MeshType>::faceGradientScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
faceGradientScheme<Type,MeshType>::faceGradientScheme
(
    const faceGradientScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
faceGradientScheme<Type,MeshType>::~faceGradientScheme()
{}

template<class Type, class MeshType>
autoPtr<faceGradientScheme<Type,MeshType>>
faceGradientScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "faceGradientSchemes", schemeName);

    word faceGradientSchemeType;
    is >> faceGradientSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(faceGradientSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face gradient scheme "
            << faceGradientSchemeType << nl << nl
            << "Valid face gradient schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceGradientScheme<Type,MeshType>>(cstrIter()(fvMsh, is));
}

}

}

}
