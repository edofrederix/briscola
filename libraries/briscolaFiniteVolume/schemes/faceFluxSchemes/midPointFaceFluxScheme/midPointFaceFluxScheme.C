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

tmp<colocatedFaceScalarField> midPointFaceFluxScheme::faceFlux
(
    const colocatedVectorField& field
)
{
    tmp<colocatedFaceScalarField> tFlux
    (
        new colocatedFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    forAllCells(Flux, i, j, k)
        Flux(i,j,k) =
            0.5
          * faceScalar
            (
                (field(i,j,k) + field(i-1,j,k)) & fan(i,j,k).left(),
                (field(i,j,k) + field(i+1,j,k)) & fan(i,j,k).right(),
                (field(i,j,k) + field(i,j-1,k)) & fan(i,j,k).bottom(),
                (field(i,j,k) + field(i,j+1,k)) & fan(i,j,k).top(),
                (field(i,j,k) + field(i,j,k-1)) & fan(i,j,k).aft(),
                (field(i,j,k) + field(i,j,k+1)) & fan(i,j,k).fore()
            );

    return tFlux;
}

tmp<staggeredFaceScalarField> midPointFaceFluxScheme::faceFlux
(
    const staggeredScalarField& field
)
{
    tmp<staggeredFaceScalarField> tFlux
    (
        new staggeredFaceScalarField
        (
            "faceFlux("+field.name()+")",
            field.fvMsh()
        )
    );

    staggeredFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    forAllDirections(Flux, d, i, j, k)
    {
        const labelVector& o = staggered::padding[d];

        const labelVector ijk(i,j,k);

        // Left/right flux always come from the x mesh direction, bottom/top
        // flux from the y mesh direction and aft/fore flux from the z mesh
        // direction. If the flux direction does not match the mesh direction,
        // we must average in the direction of mesh padding

        Flux(d,ijk) =
            0.5*fa(d,ijk)
          * faceScalar
            (
              - field(0,ijk)       - field(0,-o+ijk),
                field(0,ijk+unitX) + field(0,-o+ijk+unitX),
              - field(1,ijk)       - field(1,-o+ijk),
                field(1,ijk+unitY) + field(1,-o+ijk+unitY),
              - field(2,ijk)       - field(2,-o+ijk),
                field(2,ijk+unitZ) + field(2,-o+ijk+unitZ)
            );
    }

    return tFlux;
}

}

}

}
