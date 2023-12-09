#include "linearFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

linearFaceFluxScheme::linearFaceFluxScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dict,fvMsh)
{}

linearFaceFluxScheme::linearFaceFluxScheme
(
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dictionary(),fvMsh)
{}

tmp<colocatedLowerFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    tmp<colocatedLowerFaceScalarField> tFlux
    (
        new colocatedLowerFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedLowerFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    const colocatedFaceScalarField& fwc =
        field.fvMsh().template metrics<colocated>().faceWeightsCenter();

    const colocatedFaceScalarField& fwn =
        field.fvMsh().template metrics<colocated>().faceWeightsNeighbor();

    forAllFaces(Flux, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Flux(ijk)[fd] =
            (
                (
                    field(ijk)*fwc(ijk)[fd*2]
                  + field(nei)*fwn(ijk)[fd*2]
                )
              & fan(ijk)[fd*2]
            );
    }

    return tFlux;
}

tmp<staggeredLowerFaceScalarField> linearFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredLowerFaceScalarField> tFlux
    (
        new staggeredLowerFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    staggeredLowerFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    const staggeredFaceScalarField& fwc =
        field.fvMsh().template metrics<staggered>().faceWeightsCenter();

    const staggeredFaceScalarField& fwn =
        field.fvMsh().template metrics<staggered>().faceWeightsNeighbor();

    forAllFaces(Flux, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[d]);

        Flux(d,ijk)[fd] =
          - fa(d,ijk)[fd*2]
          * (
                field(fd,ijk)*fwc(fd,ijk)[d*2]
              + field(fd,nei)*fwn(fd,ijk)[d*2]
            );
    }

    return tFlux;
}

}

}

}
