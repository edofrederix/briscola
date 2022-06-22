#include "reconstructionScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
reconstructionScheme<Type>::reconstructionScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type>
reconstructionScheme<Type>::reconstructionScheme
(
    const reconstructionScheme<Type>& s
)
:
    scheme(s)
{}

template<class Type>
reconstructionScheme<Type>::~reconstructionScheme()
{}

template<class Type>
autoPtr<reconstructionScheme<Type>>
reconstructionScheme<Type>::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("reconstructionSchemes").subDict(name)
    );

    const word reconstructionSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reconstructionSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reconstruction scheme "
            << reconstructionSchemeType << nl << nl
            << "Valid reconstruction schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<reconstructionScheme<Type>>
    (
        cstrIter()(dict, fvMsh)
    );
}

}

}

}
