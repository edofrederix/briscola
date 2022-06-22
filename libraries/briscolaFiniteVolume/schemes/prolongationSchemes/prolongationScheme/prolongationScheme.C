#include "prolongationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
prolongationScheme<Type,MeshType>::prolongationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
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
autoPtr<prolongationScheme<Type,MeshType>> prolongationScheme<Type,MeshType>::New
(
    const word prolongationSchemeType,
    const fvMesh& fvMsh
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(prolongationSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown prolongation scheme "
            << prolongationSchemeType << nl << nl
            << "Valid prolongation schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<prolongationScheme<Type,MeshType>>
    (
        cstrIter()(dictionary(), fvMsh)
    );
}

}

}

}
