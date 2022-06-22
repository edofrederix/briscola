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
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
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
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("faceGradientSchemes").subDict(name)
    );

    const word faceGradientSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(faceGradientSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face gradient scheme "
            << faceGradientSchemeType << nl << nl
            << "Valid face gradient schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceGradientScheme<Type,MeshType>>(cstrIter()(dict, fvMsh));
}

}

}

}
