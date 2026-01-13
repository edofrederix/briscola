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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "faceFluxSchemes", schemeName);

    word faceFluxSchemeType;
    is >> faceFluxSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(faceFluxSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown face flux scheme "
            << faceFluxSchemeType << nl << nl
            << "Valid face flux schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<faceFluxScheme>(cstrIter()(fvMsh, is));
}

tmp<colocatedScalarFaceField> coloFaceFlux
(
    const staggeredScalarField& field
)
{
    tmp<colocatedScalarFaceField> tphi =
        colocatedScalarFaceField::New
        (
            "coloFaceFlux("+field.name()+")",
            field.fvMsh()
        );

    colocatedScalarFaceField& phi = tphi.ref();

    const colocatedScalarFaceField& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllFaces(phi, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        phi[fd](ijk) = -field(fd,ijk)*fa[fd](ijk);
    }

    return tphi;
}

tmp<colocatedScalarFaceField> coloFaceFlux
(
    const tmp<staggeredScalarField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<colocatedScalarFaceField> tColoFaceFlux
    (
        coloFaceFlux(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tColoFaceFlux;
}

}

}

}
