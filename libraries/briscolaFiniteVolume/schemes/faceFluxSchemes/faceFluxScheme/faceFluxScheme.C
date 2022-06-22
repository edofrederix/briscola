#include "faceFluxScheme.H"
#include "interpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

faceFluxScheme::faceFluxScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

faceFluxScheme::faceFluxScheme
(
    const faceFluxScheme& s
)
:
    scheme(s)
{}

faceFluxScheme::~faceFluxScheme()
{}

autoPtr<faceFluxScheme> faceFluxScheme::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("faceFluxSchemes").subDict(name)
    );

    const word faceFluxSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(faceFluxSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face flux scheme "
            << faceFluxSchemeType << nl << nl
            << "Valid face flux schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceFluxScheme>(cstrIter()(dict, fvMsh));
}

}

}

}
