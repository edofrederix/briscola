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
                0.5
              * faceScalar
                (
                    (f(i,j,k) + f(i-1,j,k)) & fan(i,j,k).left(),
                    (f(i,j,k) + f(i+1,j,k)) & fan(i,j,k).right(),
                    (f(i,j,k) + f(i,j-1,k)) & fan(i,j,k).bottom(),
                    (f(i,j,k) + f(i,j+1,k)) & fan(i,j,k).top(),
                    (f(i,j,k) + f(i,j,k-1)) & fan(i,j,k).aft(),
                    (f(i,j,k) + f(i,j,k+1)) & fan(i,j,k).fore()
                );
        }
    }

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

    forAll(field, l)
    {
        const staggeredScalarDirection& f0 = field[l][0];
        const staggeredScalarDirection& f1 = field[l][1];
        const staggeredScalarDirection& f2 = field[l][2];

        forAll(field[l], d)
        {
            staggeredFaceScalarDirection& F = Flux[l][d];

            const labelVector padding(staggered::padding[d]);

            const labelVector ox(d == 0 ? unitX : padding);
            const labelVector oy(d == 1 ? unitY : padding);
            const labelVector oz(d == 2 ? unitZ : padding);

            const staggeredFaceScalarDirection& fa =
                this->fvMsh().metrics<staggered>().faceAreas()[l][d];

            // Left/right flux always come from the x mesh direction, bottom/top
            // flux from the y mesh direction and aft/fore flux from the z mesh
            // direction. If the flux direction does not match the mesh
            // direction, we must average in the direction of mesh padding

            F = Zero;

            forAllCells(F, i, j, k)
            {
                const labelVector ijk(i,j,k);

                const labelVector ijkr(ijk+unitX);
                const labelVector ijkt(ijk+unitY);
                const labelVector ijkf(ijk+unitZ);

                F(ijk) =
                    0.5*fa(ijk)
                  * faceScalar
                    (
                      - f0(ijk)  - f0(ijk  - ox),
                        f0(ijkr) + f0(ijkr - ox),
                      - f1(ijk)  - f1(ijk  - oy),
                        f1(ijkt) + f1(ijkt - oy),
                      - f2(ijk)  - f2(ijk  - oz),
                        f2(ijkf) + f2(ijkf - oz)
                    );
            }
        }
    }

    return tFlux;
}

}

}

}
