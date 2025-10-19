#include "faceRestrictionScheme.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
faceRestrictionScheme<Type,MeshType>::faceRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
faceRestrictionScheme<Type,MeshType>::faceRestrictionScheme
(
    const faceRestrictionScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
faceRestrictionScheme<Type,MeshType>::~faceRestrictionScheme()
{}

template<class Type, class MeshType>
autoPtr<faceRestrictionScheme<Type,MeshType>>
faceRestrictionScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word fieldName
)
{
    return faceRestrictionScheme<Type,MeshType>::NewType
    (
        fvMsh,
        fvMsh.restrictionDict().lookupOrDefault<word>
        (
            fieldName,
            faceRestrictionScheme<Type,MeshType>::defaultScheme
        )
    );
}

template<class Type, class MeshType>
autoPtr<faceRestrictionScheme<Type,MeshType>>
faceRestrictionScheme<Type,MeshType>::NewType
(
    const fvMesh& fvMsh,
    const word faceRestrictionSchemeType
)
{
    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(faceRestrictionSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face restriction scheme "
            << faceRestrictionSchemeType << nl << nl
            << "Valid face restriction schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceRestrictionScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, scheme::nullStream)
    );
}

template<class Type, class MeshType>
void faceRestrictionScheme<Type,MeshType>::restrict
(
    faceField<Type,MeshType>& res,
    const label l
)
{
    for (label d = 0; d < MeshType::numberOfDirections; d++)
        restrict(res, l, d);
}

template<class Type, class MeshType>
void faceRestrictionScheme<Type,MeshType>::restrict
(
    faceField<Type,MeshType>& res
)
{
    for (label l = 0; l < res[0].size()-1; l++)
        restrict(res, l);
}

}

}

}
