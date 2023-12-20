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

tmp<colocatedLowerFaceScalarField> coloFaceFlux
(
    const staggeredScalarField& field
)
{
{
    tmp<colocatedLowerFaceScalarField> tphi
    (
        new colocatedLowerFaceScalarField
        (
            "coloFaceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedLowerFaceScalarField& phi = tphi.ref();

    phi = Zero;

    const colocatedFaceScalarField& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllFaces(phi, fd, i, j, k)
        phi(i,j,k)[fd] = -field(fd,i,j,k)*fa(i,j,k)[fd*2];

    return tphi;
}
}

}

}

}
