#include "upwindFaceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

upwindFaceFluxScheme::upwindFaceFluxScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dict,fvMsh)
{}

upwindFaceFluxScheme::upwindFaceFluxScheme
(
    const fvMesh& fvMsh
)
:
    faceFluxScheme(dictionary(),fvMsh)
{}

tmp<colocatedFaceScalarField> upwindFaceFluxScheme::faceFlux
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

    forAll(field, l)
    {
        colocatedFaceScalarDirection& F = Flux[l][0];

        const colocatedVectorDirection& f = field[l][0];

        const colocatedFaceVectorDirection& fan =
            this->fvMsh().metrics<colocated>().faceAreaNormals()[l][0];

        F = Zero;

        forAllCells(F, i, j, k)
        {
            F(i,j,k) =
                faceScalar
                (
                    ((f(i,j,k).x() > 0) ? f(i-1,j,k) : f(i,j,k)) & fan(i,j,k).left(),
                    ((f(i,j,k).x() > 0) ? f(i,j,k) : f(i+1,j,k)) & fan(i,j,k).right(),
                    ((f(i,j,k).y() > 0) ? f(i,j-1,k) : f(i,j,k)) & fan(i,j,k).bottom(),
                    ((f(i,j,k).y() > 0) ? f(i,j,k) : f(i,j+1,k)) & fan(i,j,k).top(),
                    ((f(i,j,k).y() > 0) ? f(i,j,k-1) : f(i,j,k)) & fan(i,j,k).aft(),
                    ((f(i,j,k).y() > 0) ? f(i,j,k) : f(i,j,k+1)) & fan(i,j,k).fore()
                );
        }
    }

    return tFlux;
}

tmp<staggeredFaceScalarField> upwindFaceFluxScheme::faceFlux
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

        if (d == 0)
        {
            Flux(d,ijk) =
                fa(d,ijk)
              * faceScalar
                (
                    ((field(0,ijk) < 0) ? - field(0,ijk) : - field(0,-o+ijk)),
                    ((field(0,ijk) > 0) ?   field(0,ijk) :   field(0, o+ijk)),
                  - 0.5 * (field(1,ijk)       + field(1,-o+ijk)),
                    0.5 * (field(1,ijk+unitY) + field(1,-o+ijk+unitY)),
                  - 0.5 * (field(2,ijk)       + field(2,-o+ijk)),
                    0.5 * (field(2,ijk+unitZ) + field(2,-o+ijk+unitZ))
                );
        }
        else if (d == 1)
        {
            Flux(d,ijk) =
                fa(d,ijk)
              * faceScalar
                (
                  - 0.5 * (field(0,ijk)       + field(0,-o+ijk)),
                    0.5 * (field(0,ijk+unitX) + field(0,-o+ijk+unitX)),
                    ((field(1,ijk) < 0) ? - field(1,ijk) : - field(1,-o+ijk)),
                    ((field(1,ijk) > 0) ?   field(1,ijk) :   field(1, o+ijk)),
                  - 0.5 * (field(2,ijk)       + field(2,-o+ijk)),
                    0.5 * (field(2,ijk+unitZ) + field(2,-o+ijk+unitZ))
                );
        }
        else if (d == 2)
        {
            Flux(d,ijk) =
                fa(d,ijk)
              * faceScalar
                (
                  - 0.5 * (field(0,ijk)       + field(0,-o+ijk)),
                    0.5 * (field(0,ijk+unitX) + field(0,-o+ijk+unitX)),
                  - 0.5 * (field(1,ijk)       + field(1,-o+ijk)),
                    0.5 * (field(1,ijk+unitY) + field(1,-o+ijk+unitY)),
                    ((field(2,ijk) < 0) ? - field(2,ijk) : - field(2,-o+ijk)),
                    ((field(2,ijk) > 0) ?   field(2,ijk) :   field(2, o+ijk))
                );
        }
    }

    return tFlux;
}

}

}

}
