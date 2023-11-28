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

tmp<colocatedFaceScalarField> coloFaceFlux
(
    const staggeredScalarField& field
)
{
{
    tmp<colocatedFaceScalarField> tphi
    (
        new colocatedFaceScalarField
        (
            "coloFaceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedFaceScalarField& phi = tphi.ref();

    phi = Zero;

    const meshField<faceScalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(phi, i, j, k)
        phi(i,j,k) =
            faceScalar
            (
              - field(0,i,  j,  k  ) * fa(i,j,k).left(),
              + field(0,i+1,j,  k  ) * fa(i,j,k).right(),
              - field(1,i,  j,  k  ) * fa(i,j,k).bottom(),
              + field(1,i,  j+1,k  ) * fa(i,j,k).top(),
              - field(2,i,  j,  k  ) * fa(i,j,k).aft(),
              + field(2,i,  j,  k+1) * fa(i,j,k).fore()
            );

    return tphi;
}
}

}

}

}
