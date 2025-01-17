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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "laplacianSchemes", schemeName);

    word laplacianSchemeType;
    is >> laplacianSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(laplacianSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown laplacian scheme "
            << laplacianSchemeType << nl << nl
            << "Valid laplacian schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<laplacianScheme<SType,Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

}

}

}
