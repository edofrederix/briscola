#include "midPointFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

midPointFaceFluxScheme::midPointFaceFluxScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dict,fvMsh)
{}

midPointFaceFluxScheme::midPointFaceFluxScheme
(
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dictionary(),fvMsh)
{}

tmp<colocatedLowerFaceScalarField> midPointFaceFluxScheme::faceFlux
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

    forAllFaces(Flux, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        Flux(ijk)[fd] = (0.5*(field(ijk) + field(nei)) & fan(ijk)[fd*2]);
    }

    return tFlux;
}

tmp<staggeredLowerFaceScalarField> midPointFaceFluxScheme::faceFlux
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

    forAllFaces(Flux, d, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[d]);

        Flux(d,ijk)[fd] =
          - 0.5*fa(d,ijk)[fd*2] * (field(fd,ijk) + field(fd,nei));
    }

    return tFlux;
}

}

}

}
