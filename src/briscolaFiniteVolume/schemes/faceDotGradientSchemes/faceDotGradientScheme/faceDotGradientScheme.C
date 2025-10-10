#include "faceDotGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
faceDotGradientScheme<Type,MeshType>::faceDotGradientScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
faceDotGradientScheme<Type,MeshType>::faceDotGradientScheme
(
    const faceDotGradientScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
faceDotGradientScheme<Type,MeshType>::~faceDotGradientScheme()
{}

template<class Type, class MeshType>
autoPtr<faceDotGradientScheme<Type,MeshType>>
faceDotGradientScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "faceDotGradientSchemes", schemeName);

    word faceDotGradientSchemeType;
    is >> faceDotGradientSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(faceDotGradientSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face gradient scheme "
            << faceDotGradientSchemeType << nl << nl
            << "Valid face gradient schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceDotGradientScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

}

}

}
