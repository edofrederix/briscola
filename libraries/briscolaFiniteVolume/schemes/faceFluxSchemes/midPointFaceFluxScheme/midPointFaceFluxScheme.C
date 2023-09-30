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
    const colocatedVectorField& f
)
{
    tmp<colocatedFaceScalarField> tFlux
    (
        new colocatedFaceScalarField
        (
            "faceFlux("+f.name()+")",
            f.fvMsh()
        )
    );

    colocatedFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const colocatedFaceVectorField& fan =
        this->fvMsh().metrics<colocated>().faceAreaNormals();

    forAllLevels(f, l, d, i, j, k)
        Flux(l,d,i,j,k) =
            0.5
          * faceScalar
            (
                (f(l,d,i,j,k) + f(l,d,i-1,j,k)) & fan(l,d,i,j,k).left(),
                (f(l,d,i,j,k) + f(l,d,i+1,j,k)) & fan(l,d,i,j,k).right(),
                (f(l,d,i,j,k) + f(l,d,i,j-1,k)) & fan(l,d,i,j,k).bottom(),
                (f(l,d,i,j,k) + f(l,d,i,j+1,k)) & fan(l,d,i,j,k).top(),
                (f(l,d,i,j,k) + f(l,d,i,j,k-1)) & fan(l,d,i,j,k).aft(),
                (f(l,d,i,j,k) + f(l,d,i,j,k+1)) & fan(l,d,i,j,k).fore()
            );

    return tFlux;
}

tmp<staggeredFaceScalarField> midPointFaceFluxScheme::faceFlux
(
    const staggeredScalarField& f
)
{
    tmp<staggeredFaceScalarField> tFlux
    (
        new staggeredFaceScalarField
        (
            "faceFlux("+f.name()+")",
            f.fvMsh()
        )
    );

    staggeredFaceScalarField& Flux = tFlux.ref();

    Flux = Zero;

    const staggeredFaceScalarField& fa =
        this->fvMsh().metrics<staggered>().faceAreas();

    forAllLevels(f, l, d, i, j, k)
    {
        const labelVector& padding = staggered::padding[d];

        const labelVector ox(d == 0 ? unitX : padding);
        const labelVector oy(d == 1 ? unitY : padding);
        const labelVector oz(d == 2 ? unitZ : padding);

        const labelVector ijk(i,j,k);

        // Left/right flux always come from the x mesh direction, bottom/top
        // flux from the y mesh direction and aft/fore flux from the z mesh
        // direction. If the flux direction does not match the mesh direction,
        // we must average in the direction of mesh padding

        Flux(l,d,ijk) =
            0.5*fa(l,d,ijk)
          * faceScalar
            (
              - f(l,0,ijk)       - f(l,0,-ox+ijk),
                f(l,0,ijk+unitX) + f(l,0,-ox+ijk+unitX),
              - f(l,1,ijk)       - f(l,1,-oy+ijk),
                f(l,1,ijk+unitY) + f(l,1,-oy+ijk+unitY),
              - f(l,2,ijk)       - f(l,2,-oz+ijk),
                f(l,2,ijk+unitZ) + f(l,2,-oz+ijk+unitZ)
            );
    }

    return tFlux;
}

}

}

}
