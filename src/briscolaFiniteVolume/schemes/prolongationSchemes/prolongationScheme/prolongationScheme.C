#include "prolongationScheme.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
prolongationScheme<Type,MeshType>::prolongationScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
prolongationScheme<Type,MeshType>::prolongationScheme
(
    const prolongationScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
prolongationScheme<Type,MeshType>::~prolongationScheme()
{}

template<class Type, class MeshType>
autoPtr<prolongationScheme<Type,MeshType>>
prolongationScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word fieldName
)
{
    return prolongationScheme<Type,MeshType>::NewType
    (
        fvMsh,
        fvMsh.prolongationDict().lookupOrDefault<word>
        (
            fieldName,
            prolongationScheme<Type,MeshType>::defaultScheme
        )
    );
}

template<class Type, class MeshType>
autoPtr<prolongationScheme<Type,MeshType>>
prolongationScheme<Type,MeshType>::NewType
(
    const fvMesh& fvMsh,
    const word prolongationSchemeType
)
{
    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(prolongationSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown prolongation scheme "
            << prolongationSchemeType << nl << nl
            << "Valid prolongation schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<prolongationScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, scheme::nullStream)
    );
}

}

}

}
